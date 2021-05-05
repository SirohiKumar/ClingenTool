from __future__ import print_function
from flask import Flask, request, render_template
from intermine.webservice import Service
service = Service("http://www.mousemine.org/mousemine/service")
import pandas as pd
from collections import OrderedDict

## this is a function that takes a gene and returns a list of MouseMine queries that can be interpreted by a separate function
## parameter is a gene symbol
def create_queries(gene):

    ## create starter queries that search for allele-based annotations
    alleleQuery = service.new_query("OntologyAnnotation")
    alleleQuery.add_constraint("subject", "Allele")
    alleleQuery.add_constraint("ontologyTerm", "DOTerm")

    ## mouse phenotypes asspciated with the gene 
    phenQuery = service.new_query("OntologyAnnotation")
    phenQuery.add_constraint("subject", "Gene")
    phenQuery.add_constraint("ontologyTerm", "MPTerm")

    ## diseases associated with the gene
    disQuery = service.new_query("OntologyAnnotation")
    disQuery.add_constraint("ontologyTerm", "DOTerm")

    queryList = [alleleQuery, phenQuery, disQuery]
    
    for query in queryList:
        ## get specific info about the publications
        query.add_view("ontologyTerm.name", "subject.symbol",
                       "evidence.publications.pubMedId", "evidence.publications.mgiId",
                       "evidence.publications.firstAuthor", "evidence.publications.title",
                       "evidence.publications.year", "evidence.publications.journal",
                       "evidence.publications.abstractText")

        ## logic (for not allele): make sure the gene is the correct one, ensure the publication is a PubMed one
        if query is not alleleQuery:
            query.add_constraint("subject.symbol", "=", gene, code = "A") 
            query.add_constraint("evidence.publications.pubMedId", "IS NOT NULL", code = "B") 
            query.set_logic("A and B")

    ## allele logic: make sure the publication is from PubMed, is the correct gene name, and 
    alleleQuery.add_constraint("evidence.publications.pubMedId", "IS NOT NULL", code = "A") 
    alleleQuery.add_constraint("subject.feature.symbol", "=", gene, code = "B")
    alleleQuery.add_constraint("evidence.baseAnnotations", "IS NULL", code = "C")
    alleleQuery.set_logic("A and B and C")

    return queryList

## this function runs a list of queries and then outputs the data into a dict
## parameters are a list of queries (from prev function) and a dictionary to add them to
def run_queries(queryList, geneInfoDict):
    
    ## determine which query has phenotype information THIS DEPENDS ON THE ORDER OF THE QUERIES
    phenQuery = queryList[1]

    ## for each query...
    for query in queryList:

        ## mark whether it has disease or phenotype information 
        isDisease = True
        if query is phenQuery:
            isDisease = False

        ## ...loop through all the publications that fulfull its conditions
        for row in query.rows():

            ## pull the PubMed ID to check if this publication is already in the dictionary
            currentPub = str(row["evidence.publications.pubMedId"])

            ## if it's not already in the dictionary, give it its own entry + add the relevant publication details
            if currentPub not in geneInfoDict:
                geneInfoDict[currentPub] = {}

                geneInfoDict[currentPub]["Author"] = str(row["evidence.publications.firstAuthor"])
                geneInfoDict[currentPub]["Title"] = str(row["evidence.publications.title"])
                geneInfoDict[currentPub]["Journal"] = str(row["evidence.publications.journal"])
                geneInfoDict[currentPub]["Year"] = str(row["evidence.publications.year"])
                geneInfoDict[currentPub]["Abstract"] = str(row["evidence.publications.abstractText"])

                ## ft empty lists for disease and phenotype information that hasn't been found yet
                geneInfoDict[currentPub]["Diseases"] = []
                geneInfoDict[currentPub]["Phenotypes"] = []

            ## pull the query's annotation: associated phenotype or disease 
            data  = str(row["ontologyTerm.name"])

            ## if the query has disease information add that to the list of disease information
            if isDisease == True:
                dataList = geneInfoDict[currentPub]["Diseases"]
                dataList.append(data)
                
                geneInfoDict[currentPub]["Diseases"] = dataList

            ## if the query has phenotype information, add it to the list of phen information
            if isDisease == False:
                dataList = geneInfoDict[currentPub]["Phenotypes"]
                dataList.append(data)
                
                geneInfoDict[currentPub]["Phenotypes"] = dataList

    ## return the updated geneInfoDict
    return geneInfoDict

## function to check if a gene is in MGI (to avoid error massages).
## parameter is a gene symbol
def check_gene(gene):
    query = service.new_query("Gene")
    query.add_view("publications.pubMedId")
    query.add_constraint("symbol", "=", gene, code = "A")

    x = 0
    for row in query.rows():
        x = x + 1
    ## if the query says there are no associated publications in MGI...
    if x == 0:
        return False
    ## if there are associated publications...
    else:
        return True

app = Flask(__name__)

@app.route('/', methods=['POST', 'GET'])
def hello_world():
    ## if the person has been here before aka they hit the submit button and need results
    if request.method == 'POST':
        gene = request.form['gene']
        phen = request.form['phen']
        dis = request.form['dis']
        rescue = request.form['rescue']
        includeAbstract = request.form['abstract']
        removeDis = []
        removePhen = []
        yesRes = []
        noRes = []

        ## all information for the user
        fields = ["Author", "Title", "Journal", "Year", "Abstract", "Diseases", "Phenotypes"]

        ## terms used to determine if a publication is rescue or not
        rescueTerms = ['rescue', 'rescued', 'rescues', 'ameliorate', 'ameliorated', 'ameliorates', 'transduction', 'transduces', 'transduce', 'complementation']

        ## get rid of extra spaces in the gene name and upper-case it (queries are NOT case-sensitive but)
        gene = gene.replace(" ", "")
        gene = gene.upper()

        ## get essential data
        geneInfoDict = {}

        queryList = create_queries(gene)
        geneInfoDict = run_queries(queryList, geneInfoDict)
        
        ## ERROR MESSAGE IF THE GENE ENTERED IS NON-EXISTENT
        if(check_gene(gene) == False):
            errorMessage = render_template('error.html', gene=gene, message="This gene does not exist in the MGI database")
            return errorMessage
        

        ## IF THERE ARE RESULTS THEN RETURN A FORMATTED DICTIONARY
        if (check_gene(gene) == True):
            for currentPub in geneInfoDict:
                
                ## REPLACE W/ STRING VERSION OF DISEASE LIST
                diseaseList = geneInfoDict[currentPub]["Diseases"]
                ##remove duplicates + sort list alphabetically
                diseaseList = list(set(diseaseList))
                diseaseList.sort()
                strDiseaseList = ""
                diseaseCount = len(diseaseList)
                ## if there are no diseases:  
                if diseaseCount == 0:
                    strDiseaseList = "No associated diseases"
                ## if there is only one disease
                if diseaseCount == 1:
                    strDiseaseList = str(diseaseList[0])
                    strDiseaseList = strDiseaseList.capitalize()
                ## if there are multiple diseases: 
                if diseaseCount > 1:
                    for disease in diseaseList:
                        if disease == diseaseList[0]:
                            strDiseaseList = strDiseaseList + disease.capitalize() + ", "
                        elif disease == diseaseList[(diseaseCount - 1)]:
                            strDiseaseList = strDiseaseList + disease
                        else:
                            strDiseaseList = strDiseaseList + disease + ", "
                ## replace the list in the dictionary
                geneInfoDict[currentPub]["Diseases"] = strDiseaseList

                ## REPLACE WITH STRING VERSION OF PHENOTYPE LIST
                phenotypeList = geneInfoDict[currentPub]["Phenotypes"]
                stringPhenotypeList = ""
                ## remove duplicate phenotypes + sort the list alphabetically
                phenotypeList = list(set(phenotypeList))
                phenotypeList.sort()
                numberOfPhen = len(phenotypeList)
                ## turn the list of phenotypes into a string version of itself
                if numberOfPhen == 0:
                    stringPhenotypeList = "No associated phenotypes"
                if numberOfPhen == 1:
                    phenotype = phenotypeList[0]
                    ## phenotype formatting to make it less confusing
                    phenotype = phenotype.replace(", complete penetrance", " (complete penetrance)")
                    phenotype = phenotype.replace(", incomplete penetrance", " (incomplete penetrance)")
                    stringPhenotypeList = stringPhenotypeList + phenotype.capitalize()
                if numberOfPhen > 1:
                    lastPhenotype = phenotypeList[numberOfPhen - 1]
                    firstPhen = phenotypeList[0]
                    for phenotype in phenotypeList:
                        phenotype = phenotype.replace(", complete penetrance", " (complete penetrance)")
                        phenotype = phenotype.replace(", incomplete penetrance", " (incomplete penetrance)")
                        if phenotype == firstPhen:
                            stringPhenotypeList = stringPhenotypeList + firstPhen.capitalize() + ", "
                        elif phenotype != lastPhenotype:
                            stringPhenotypeList = stringPhenotypeList + phenotype + ", "
                        if phenotype == lastPhenotype:
                            stringPhenotypeList = stringPhenotypeList + phenotype
        
                geneInfoDict[currentPub]["Phenotypes"] = stringPhenotypeList

                ## collect publications without associated DISEASES IF INDICATED
                if dis == "yes":
                    if geneInfoDict[currentPub]["Diseases"] == "No associated diseases":
                        removeDis.append(currentPub)

                ## collect publications without associated PHENOTYPES IF INDICATED
                if phen == "yes":
                    if geneInfoDict[currentPub]["Phenotypes"] == "No associated phenotypes":
                        removePhen.append(currentPub)

                ## FILTER OUT RESCUE PAPERS IF INDICATED (tentative) IF INDICATED
                if rescue == "yes":
                    fields.append("Is this a rescue paper?")
                    abstract = geneInfoDict[currentPub]["Abstract"]
                    abstractList = abstract.split()
                    strAbs = ""
                    ## check if any of the words match the rescue words
                    for word in abstractList:
                        newWord = word.lower()
                        ## remove potentially obstructive punctuation
                        newWord = newWord.replace(",", "")
                        newWord = newWord.replace(".", "")
                        newWord = newWord.replace("(", "")
                        newWord = newWord.replace(")", "")
                        ## check if the word is in the list of rescue terms
                        if newWord in rescueTerms:
                            strAbs = strAbs + " " + "stBold"
                            strAbs = strAbs + word
                            strAbs = strAbs + "clBold"
                            yesRes.append(currentPub)
                        else:
                            strAbs = strAbs + " " + word
                    ## if the publication isn't a rescue paper, add it to a list
                    if currentPub not in yesRes:
                        noRes.append(currentPub)

                    geneInfoDict[currentPub]["Abstract"] = strAbs

                ## INDICATE IF A PUBLICATION IS RESCUE OR NOT
                for pmid3 in yesRes:
                    geneInfoDict[pmid3]["Is this a rescue paper?"] = "Yes"
                for pmid4 in noRes:
                    geneInfoDict[pmid4]["Is this a rescue paper?"] = "No"

                ## IF THERE IS NO AVAILABLE ABSTRACT TELL THE USER
                abstract = geneInfoDict[currentPub]["Abstract"]
                if abstract == "None":
                    abstract = "No abstract available for this publication"

                geneInfoDict[currentPub]["Abstract"] = abstract

                ## MARK HTML TAGS THAT WERE ALREADY IN THE TITLES + ABSTRACTS WHEN PULLED FROM MGI
                    # the to.html() function treats each cell as a whole string and the tags need to be 
                    # added after the function has been used on the dataframe
                title = geneInfoDict[currentPub]["Title"]

                title = title.replace('<i>', 'stItal')
                title = title.replace('</i>', 'enItal')
                title = title.replace('<sup>', 'stSup')
                title = title.replace('</sup>', 'enSup')

                geneInfoDict[currentPub]["Title"] = title

                ## get rid of the ABSTRACT IF INDICATED + IF RESCUE IS A NO
                ## ANYTHING HAVING TO DO WITH ABSTRACTS SHOULD BE CODED BEFORE HERE
                if includeAbstract == "no" and currentPub not in yesRes:
                    geneInfoDict[currentPub]["Abstract"] = " "

                ## ENSURE ORDER OF DICT IS MAINTAINED
                currentDict = OrderedDict(geneInfoDict[currentPub])
                for header in fields:
                    currentDict.move_to_end(header)
                geneInfoDict[currentPub] = dict(currentDict)

            ## get rid of unwanted publications w/ different countervariables (to avoid confusion)
            for pmid1 in removeDis:
                geneInfoDict.pop(pmid1)
            for pmid2 in removePhen:
                geneInfoDict.pop(pmid2)

            pubCount = len(geneInfoDict)

            ## turn into a dataframe w/ pandas
            geneInfoDF = pd.DataFrame(data=geneInfoDict)
            ## reverse columns and rows
            geneInfoDF = geneInfoDF.T

            ## PULL RESULTS PAGE TEMPLE AND MODIFY BELOW
            geneInfoHTML = render_template('results.html')

            geneInfoHTML = geneInfoHTML.replace("@@@@", geneInfoDF.to_html(classes='mystyle'))

            # FORMAT THE TABLE

            # replace the PMIDs in the dataframe with html so they read as cells instead of headers
            for pmid in geneInfoDict:
                pmid = str(pmid)
                originalHeader = "<th>"+pmid+"</th>"
                link = "'https://pubmed.ncbi.nlm.nih.gov/"+pmid+"/'"
                newCell = "<td><a href="+link+"target='_blank'>"+pmid+"</a></td>"
                geneInfoHTML = geneInfoHTML.replace(originalHeader, newCell)

            ## add in the correct table id for the script
            geneInfoHTML = geneInfoHTML.replace('<table border="1" class="dataframe mystyle">', '<table id="resultsTable" border="1" class="dataframe mystyle">')

            ## replace headers with "clickable" headers (order is important)
            geneInfoHTML = geneInfoHTML.replace('<th></th>', '<th onclick="sortTable(0)" title="Click to sort">PMID</th>')
            geneInfoHTML = geneInfoHTML.replace('<th>Author</th>', '<th onclick="sortTable(1)" title="Click to sort">Author</th>')
            geneInfoHTML = geneInfoHTML.replace('<th>Title</th>', '<th onclick="sortTable(2)" title="Click to sort">Title</th>')
            geneInfoHTML = geneInfoHTML.replace('<th>Journal</th>', '<th onclick="sortTable(3)" title="Click to sort">Journal</th>')
            geneInfoHTML = geneInfoHTML.replace('<th>Year</th>', '<th onclick="sortTable(4)" title="Click to sort">Year</th>')
            if includeAbstract == "yes" or rescue == "yes":
                geneInfoHTML = geneInfoHTML.replace('<th>Abstract</th>', '<th onclick="sortTable(5)" title="Click to sort">Abstract</th>')
            geneInfoHTML = geneInfoHTML.replace('<th>Diseases</th>', '<th onclick="sortTable(6)" title="Click to sort">Diseases</th>')
            geneInfoHTML = geneInfoHTML.replace('<th>Phenotypes</th>', '<th onclick="sortTable(7)" title="Click to sort">Phenotypes</th>')

            ## if rescue papers also click on header
            if rescue == "yes":
                newHeader = '<th onclick="sortTable(8)" title="Click to sort (twice for desc)">Is this a rescue paper?</th>'
                geneInfoHTML = geneInfoHTML.replace('<th>Is this a rescue paper?</th>', newHeader)
                rescueString = '''<br>

                <h1>Papers are marked as rescue if their abstracts include key terms (list below). These terms are marked in the abstracts of these rescue publications. </h1>
                
                <label>Terms:</label>
                <li>Rescue</li>
                <li>Ameliorate</li>
                <li>Transduce</li>
                <li>Complementation</li>
                <br>
                '''
                geneInfoHTML = geneInfoHTML.replace('[rescue]', rescueString)

            if rescue == "no":
                geneInfoHTML = geneInfoHTML.replace('[rescue]', '')


            ## REPLACE HTML INDICATORS INSERTED ABOVE WITH FUNCTIONING HTML TAGS
            geneInfoHTML = geneInfoHTML.replace('stBold', '<b><mark>')
            geneInfoHTML = geneInfoHTML.replace('clBold', '</mark></b>')
            geneInfoHTML = geneInfoHTML.replace('stItal', '<i>')
            geneInfoHTML = geneInfoHTML.replace('enItal', '</i>')
            geneInfoHTML = geneInfoHTML.replace('stSup', '<sup>')
            geneInfoHTML = geneInfoHTML.replace('enSup', '</sup>')

            ## replace other placeholders with info abt the dict
            geneInfoHTML = geneInfoHTML.replace("****", gene)
            geneInfoHTML = geneInfoHTML.replace("$$$$", str(pubCount))

            ## IF THERE ARE 0 PUBLICATIONS LEFT RETURN THE (MODIFIED) ERROR MESSAGE: 
            if pubCount == 0:
                message = "There are currently 0 publications that show an association between this gene and any diseases or phenotypes"
                errorMessage = render_template('error.html', gene=gene, message=message)
                return errorMessage

            if pubCount == 1:
                geneInfoHTML = geneInfoHTML.replace("You are viewing all 1 publications", "You are viewing the 1 publication")
                return geneInfoHTML

            ## IF THERE ARE RESULTS, PRINT THE TABLE
            else:
                return geneInfoHTML

    ## if the person hasn't visited the page before yet, ask them to enter a gene name aka
    ## the original page people see when they visit the website
    else:
        ## pull form from templates folder
        return render_template('form.html')







