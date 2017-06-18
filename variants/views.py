from __future__ import unicode_literals
from django.shortcuts import get_object_or_404, redirect
from .models import Annotation,Population,Variant
from django.shortcuts import render,render_to_response
from django.conf import settings
from django.db.models import Q
import json

import re
import time
import os
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt
from django.template import  loader
import sys
from django.db import connection


#input : String from the user
#

##new
#Query = {"query":["Gene_Symbol","NRIP1","AND","Allele","AAGGTCCATTGA"]} # :'D

##Query = {"query":["Gene_Symbol","NRIP1", "OR", "Position", "9825793"]}

### test
Q_list = []

Query1 = {"query":["Allele_Count","29", "AND", "Filter", "PASS", "AND", "Position", "37858500"]}
Query2 = {"query":["Filter","Pass", "OR", "Allele", "TGG", "OR", "AC_Adjusted", "0"]}
Query3 = {"query":["Reference","G", "AND", "Gene_Symbol", "MIR3687"]}
Query4 = {"query":["AC_Adjusted","0", "AND", "Feature", "Transcript"]}
Query5 = {"query":["Impact","MODERATE", "AND", "Filter", "InbreedingCoeff_Filter"]}
Query6 = {"query":["Reference","TG", "OR", "Consequence", "upstream_gene_variant", "OR", "AC_Adjusted", "1987"]}
Query7= {"query":["Filter","InbreedingCoeff_Filter", "AND", "Chromosome_Number", "21", "AND", "AN_Adjusted", "0"]}
Query8 = {"query":["Alternate","AGG", "AND", "Gene_Symbol", "TEKT4P2"]}
Query9 = {"query":["Filter","Pass", "AND", "Position", "9907203", "OR", "AN_Adjusted", "44"]}
Query10 = {"query":["Filter","Pass", "AND", "Impact", "HIGH"]}
Query11 = {"query":["HGVS_Protein","ENST00000609713.1:c.999T>C(p.%3D)", "OR", "Reference", "TTCCAGACTGGGTACA"]}
Query12 = {"query":["Position","9825838", "AND", "Consequence", "synonymous_variant"]}
Query13 = {"query":["AC_Adjusted","0", "AND", "Gene_Symbol", "TEKT4P2"]}
Query14 = {"query":["AN_Adjusted","0", "AND", "Chromosome_Number", "21", "AND", "Impact", "HIGH"]}
Query15 = {"query":["Filter","Pass", "AND", "Impact", "LOW", "AND", "Position", "9907224"]}
Query16= {"query":["Position","9907237", "OR", "Gene", "ENSG00000273254"]}
Query17 = {"query":["PolyPhen","benign(0.005)", "OR", "Reference", "TTCC", "AND", "Allele", "T"]}
Query18 = {"query":["Filter","AC_Adj0_Filter", "AND", "Consequence", "start_lost", "OR", "Gene", "ENSG00000160255"]}
Query19 = {"query":["Gene_Symbol","SON", "AND", "Consequence", "synonymous_variant"]}
Query20 = {"query":["Gene_Symbol","BRWD1", "AND", "Impact", "HIGH", "AND", "Consequence", "synonymous_variant"]}
Query21 = {"query":["SIFT","tolerated_low_confidence(0.13)", "AND", "PolyPhen", "possibly_damaging(0.702)"]}
Query22 = {"query":["Position","9907224", "OR", "HGVS_Protein", "ENST00000609713.1:c.999T>C(p.%3D)"]}

Q_list.append(Query1)
Q_list.append(Query2)
Q_list.append(Query3)
Q_list.append(Query4)
Q_list.append(Query5)
Q_list.append(Query6)
Q_list.append(Query7)
Q_list.append(Query8)
Q_list.append(Query9)
Q_list.append(Query10)
Q_list.append(Query11)
Q_list.append(Query12)
Q_list.append(Query13)
Q_list.append(Query14)
Q_list.append(Query15)
Q_list.append(Query16)
Q_list.append(Query17)
Q_list.append(Query18)
Q_list.append(Query19)
Q_list.append(Query20)
Q_list.append(Query21)
Q_list.append(Query22)


###


#Query_list_input = {"Query_list": [{"column_name": "Gene_Symbol", "value": "NRIP1", "operation": "AND"}, {"column_name": "Allele", "value": "AAGGTCCATT", "operation": "null"}]}


## JUST IN CASE ######################################
#
# Variant_attributes = ["Chromosome_Number", "Position", "Var_ID", "Reference", "Alternate", "Quality", "Filter"]
# Population_attributes = ["Variant",  "Allele_Count",  "AC_African_American",  "AC_American",  "AC_Adjusted",  "AC_East_Asian",  "AC_Finnish",  "AC_Hemizygous",  "AC_Heterozygous",  "AC_Homozygous",  "AC_Non_Finnish_European",  "AC_Other",  "AC_South_Asian",  "Allele_Frequency",  "Allele_Number",  "AN_African_American",  "AN_American",  "AN_Adjusted",  "AN_East_Asian",  "AN_Finnish",  "AN_Non_Finnish",  "AN_Other",  "AN_South_Asian",  "Read_Depth",  "Genotype_Quality_MEAN",  "GQ_Standard_Deviation",  "Hemi_African_American",  "Hemi_American",  "Hemi_East_Asian",  "Hemi_Finnish",  "Hemi_Non_Finnish_European",  "Hemi_Other",  "Hemi_South_Asian",  "Het_African_American",  "Het_American",  "Het_East_Asian",  "Het_Finnish",  "Het_Non_Finnish_European",  "Het_Other",  "Het_South_Asian",  "Hom_African_American",  "Hom_American",  "Hom_East_Asian",  "Hom_Finnish",  "Hom_Non_Finnish_European_Homozygous",  "Hom_Other",  "Hom_South_Asian",  "InbreedingCoef",  "DP_HIST",  "GQ_HIST"]
# Annotation_attributes = ["Variant",  "Population",  "Allele",  "Consequence",  "Impact",  "Gene_Symbol",  "Gene",  "Feature",  "Feature_type",  "Biotype",  "Exon",  "Intron",  "HGVS_DNA",  "HGVS_Protein",  "cDNA_position",  "CDS_position",  "Protein_position",  "Amino_acids",  "Codons",  "Existing_variation",  "ALLELE_NUM",  "DISTANCE",  "STRAND",  "VARIANT_CLASS",  "MINIMISED",  "SYMBOL_SOURCE",  "HGNC_ID",  "CANONICAL",  "TSL",  "CCDS",  "ENSP",  "SWISSPROT",  "TREMBL",  "UNIPARC",  "SIFT",  "PolyPhen",  "DOMAINS",  "HGVS_OFFSET",  "GMAF",  "AFR_MAF",  "AMR_MAF",  "ASN_MAF",  "EAS_MAF",  "EUR_MAF",  "SAS_MAF",  "AA_MAF",  "EA_MAF",  "CLIN_SIG",  "SOMATIC",  "PHENO",  "PUBMED",  "MOTIF_NAME",  "MOTIF_POS",  "HIGH_INF_POS",  "MOTIF_SCORE_CHANGE",  "LoF_info",  "LoF_flags",  "LoF_filter",  "LoF",  "context",  "ancestral"]
#
#
# def attribute_exists(in_attribute):
#
#     if hasattr(Variant, in_attribute):
#         return "exists in Variant"
#     elif hasattr(Population, in_attribute):
#         return "exists in Population"
#     elif hasattr(Annotation, in_attribute):
#         return "exists in Annotation"
#
#     return "Does not exist"
#######################################################



def Complex_query(index, Query_list):

    join_string = ""

    if hasattr(Variant, Query_list[index]['column_name']):
        join_string = "Variant__"
    elif hasattr(Annotation, Query_list[index]['column_name']):
        join_string = ""
    elif hasattr(Population, Query_list[index]['column_name']):
        join_string = "Population__"
    else :
        print"Attribute does not exist in any table"
        #sys.exit(1)


    if (Query_list[index]['operation'] == "null"):
        return Q(**{join_string + Query_list[index]['column_name']: Query_list[index]['value']})

    if (Query_list[index]['operation'] == "AND"):
        return Q(**{join_string + Query_list[index]['column_name']: Query_list[index]['value']}) & Complex_query(index + 1, Query_list)

    if (Query_list[index]['operation'] == "OR"):
        return Q(**{join_string + Query_list[index]['column_name']: Query_list[index]['value']}) | Complex_query(index + 1, Query_list)


##parser_start

def Parse(Query_input_str):

    query_list = Query_input_str['query']

    #parsed_str = "{\"Query_list\": ["
    parsed_str = "["

    i = 0;
    while True:
        parsed_str += "{"
        parsed_str += "\"column_name\": \""
        parsed_str += query_list[i]
        i += 1
        parsed_str += "\", \"value\": \""
        parsed_str += query_list[i]
        i += 1
        parsed_str += "\", \"operation\": \""

        if i >= len(query_list):
            parsed_str += "null"
            parsed_str += "\"}"
            break

        parsed_str += query_list[i]
        i += 1
        parsed_str += "\"}, "

    #parsed_str += "]}"
    parsed_str += "]"

    #Query_list_input = json.loads(parsed_str)
    #Query_list = Query_list_input['Query_list']

    Query_list = json.loads(parsed_str)

    return Query_list

##parser_end



##newend
def get_type(line):

    # Defining patterns to be identified using regular expressions
    Gene = re.compile('^G=.+')
    Transcript = re.compile('^T=.+')
    Variant = re.compile('^V=.+')
    Region = re.compile('^R=.+')
    AlleleFrequency = re.compile('^AF=.+')
    Population = re.compile('^POP=.+')

    ComplexQuery = re.compile('^CQ=.+')

    ComplexQuerydebug = re.compile('^debug=.+')

    # Empty Pattern
    Empty = re.compile('^.*=$')

    # Matching Input String to a Specific pattern from the pre-defined one.
    if Gene.match(line):
        return "Gene"
    elif Variant.match(line):
        return "Variant"
    elif Transcript.match(line):
        return "Transcript"
    elif Region.match(line):
        return "Region"
    elif AlleleFrequency.match(line):
        return "AlleleFrequency"
    elif Population.match(line):
        return "Population"
    elif ComplexQuery.match(line):
        return "ComplexQuery"
    elif ComplexQuerydebug.match(line):
        return "ComplexQueryDebug"
    elif Empty.match(line):
        return "Empty"

    #If No Match Found
    else:
        print 'NOT'


##############################[ROUTER:START]############################
def Route(request):
    if 'Query' in request.POST:
        #
        # return redirect('/variants/Gene_view/{}'.format(request.POST['Query']))
        #
        # request.session['_old_post'] = request.POST
        # return HttpResponseRedirect('Gene_view')
        #
        # Get Input Posted from HTML Template

        print request.POST['Query']
        type = get_type(request.POST['Query'])

        # Compare it's type to Given Types Then Redirect to the Right View
        if (type == "Gene"):
            return Gene_view(request)
        elif (type == "Variant"):
            return variant_view(request)
            #return variant_view(request)
        elif (type == "Transcript"):
            return Transcript_view(request)
        elif (type == "Region"):
            return Region_view(request)
        elif (type == "ComplexQuery"):
            return ComplexQuery_view(request)
        elif (type == "ComplexQueryDebug"):
            print "heloo"
            return ComplexQuerydebug_view(request)
        else:
            return error_view(request)
    else:
        return
##############################[ROUTER:END]##############################


#HOME_VIEW

def Home_view(request):

    return render(request,'Home_view.html')


def Advanced_Search_view(request):
    return render(request, 'Advanced_Search.html')


def error_view(request):
    return render(request, 'error.html')


def variant_view(request):
    start_time=time.time();
    genes=[]
    transcripts=[]

    if 'Query' in request.POST:
        vr_id=request.POST['Query'].split('=', 1)[-1]
        VARIANT=str(vr_id)
        try:
            # for variant input 1:13723 C/G
            VARIANT = VARIANT.replace(':', '-')
            VARIANT = VARIANT.replace(' ', '-')
            VARIANT = VARIANT.replace('/', '-')

            variant_list=VARIANT.split("-")

            # get the query variant by Chromosome number,Position & Alternate
            variant = get_object_or_404(Variant,Chromosome_Number=variant_list[0],Position=variant_list[1],Alternate=variant_list[3])
            vr_id=variant.id

            # get population record of this variant
            varpop = get_object_or_404(Population,Variant_id=vr_id)

            # calculate allele frequency for each population and Total
            allele_freq=calc_allele_freq(varpop)

            # get annotation records of the query variant
            var_anno = Annotation.objects.filter(Variant_id=vr_id)

            # for displaying annotation in specific format ,for more details Go the function declaration
            list_of_annotations = Display_annotations(var_anno)

            # Counting Transcripts and Genes in the query variant
            for var in var_anno:
                transcripts.append(var.Feature_type)
                genes.append(var.Gene_Symbol)

            #for eliminating empty strings of genes
            genes = filter(None, genes)

            trans_num = len(set(transcripts))
            gene_num = len(set(genes))

            query_time = format(float(time.time() - start_time), '.4g')

            # Saving time of query to a file
            target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
            target.write("Variant=")
            target.write(VARIANT)
            target.write("  ")
            target.write(query_time)
            target.write("\n")
            target.close()
            return render(request, 'Variant_view.html',
                          {'variant': variant,'varpop' :varpop,'allele_freq':allele_freq,'trans_num':trans_num,
                           'gene_num':gene_num ,'list_of_annotations':list_of_annotations,'query_time':query_time})
        except:
            return error_view(request)
    else:
        return



def Gene_view(request):
    start_time = time.time();
    ID_list = []
    Transcript_list = []
    if 'Query' in request.POST:
        gene_symbol = request.POST['Query'].split('=', 1)[-1]
        test_var = str(gene_symbol)

        # Check if Request is Gene symbol or Gene_code to Determine Filter Type
        try:
            if test_var[:4] == 'ENSG':
                Annotations = Annotation.objects.filter(Gene=gene_symbol)
                gene_symbol = Annotations[0].Gene_Symbol
            else:
                Annotations = Annotation.objects.filter(Gene_Symbol=gene_symbol)

                # Annotations = Annotation.objects.filter(Gene_Symbol=gene_symbol)
            gene_ID = Annotations[0].Gene

            # Retrieves Unique Foreign Keys & Transcripts of each GENE
            for rec in Annotations:
                ID_list.append(rec.Variant_id)
                Transcript_list.append(rec.Feature_type)

            Transcripts = set(Transcript_list)
            var_id_list = set(ID_list)

            Variants, filtered_variants, Genes_ID = Display_variants(var_id_list)
            query_time = format(float(time.time() - start_time), '.4g')

            # Saving time of query to a file
            target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
            target.write("Gene=")
            target.write(test_var)
            target.write("  ")
            target.write(query_time)
            target.write("\n")
            target.close()
            return render(request, 'Gene_view.html',
                          {'gene_ID': gene_ID, 'Variants_num': len(var_id_list), 'gene_symbol': gene_symbol,
                           'Variants': Variants, 'filtered_variants': filtered_variants, 'Transcripts': Transcripts,
                           'query_time': query_time})
        except:
            return error_view(request)
    else:
        return request



def Transcript_view(request):
    start_time = time.time();
    #global Gene_symbol
    ID_list = []
    Transcript_list=[]
    if'Query'in request.POST:
        transcript = request.POST['Query'].split('=', 1)[-1]
        try:
            # Filter each Annotation by Transcript
            Annotations = Annotation.objects.filter(Feature_type=transcript)

            # Retrieves Gene symbol ex:MIR3687 of Transcript if exists
            Gene_symbol = Annotations[0].Gene_Symbol

            # Retrieves other Transcripts in same GENE & Unique Foreign Keys of all Annotations

            for ann in Annotations:
                ID_list.append(ann.Variant_id)
                Transcript_list.append(ann.Feature_type)

            Transcripts = set(Transcript_list)
            var_id_list = set(ID_list)

            Variants , filtered_variants, Genes = Display_variants(var_id_list)
            query_time = format(float(time.time() - start_time), '.4g')

            # Saving time of query to a file
            target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
            target.write("Transcript=")
            target.write(transcript)
            target.write("  ")
            target.write(query_time)
            target.write("\n")
            target.close()
            return render(request, 'Transcript_view.html', {'Variants': Variants, 'Gene_symbol': Gene_symbol, 'transcript': transcript, 'Variants_num': len(var_id_list),
                                                            'filtered_variants': filtered_variants, 'Transcripts': Transcripts,'query_time':query_time })
        except:
            return error_view(request)



def Region_view(request):
    start_time = time.time();
    ID_list = []

    if 'Query' in request.POST:
        region = request.POST['Query'].split('=', 1)[-1]
        Test_var = str(region)
        try:
            #Split the query region by input"-"character
            Region_split = Test_var.split("-")
            #Display region in format chrom/position1/position2
            Region_Info = str(Region_split[0])
            Region_Info += " / "
            Region_Info += str(Region_split[1])
            Region_Info += " / "
            Region_Info += str(Region_split[2])

            #get all variants  between position1 and position2
            Variant_list = Variant.objects.filter(Chromosome_Number=Region_split[0], Position__range=(Region_split[1], Region_split[2]))
            for var in Variant_list:
                ID_list.append(var.id)
            # retrieving list of all variants with it's info , count of filtered varianats and list of Genes_ID
            Variants, filtered_variants, Genes_ID = Display_variants(ID_list)

            query_time = format(float(time.time() - start_time), '.4g')

            # Saving time of query to a file
            target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
            target.write("Region=")
            target.write(region)
            target.write("  ")
            target.write(query_time)
            target.write("\n")
            target.close()
            return render(request, 'Region_view.html', {'Genes_ID': Genes_ID, 'Region_Info': Region_Info, 'Variants': Variants,'query_time':query_time})

        except:
            return error_view(request)



def Display_annotations(object):
    # each item in list_of_annotations has [1 annotation,1/more list of genes_transcripts]
    # genes_transcripts = [gene,[transcripts]]
    list_of_annotations = []
    transcripts = []
    genes_transcripts=[]
    # Append all annotations as a item_list in list_of_annotations
    for item in object:
        temp = []
        temp.append(item.Consequence)
        list_of_annotations.append(temp)

    #Unique_annotations
    list_of_annotations=set(map(tuple,list_of_annotations))
    list_of_annotations=map(list,list_of_annotations)

    #Append all genes have the same annotation in genes_transcripts list as item_list
    for item in list_of_annotations:
        for i in object:
            if item[0] == i.Consequence:
                gene = []
                gene.append(i.Gene_Symbol)
                genes_transcripts.append(gene)
        #Unique genes
        genes_transcripts=set(map(tuple,genes_transcripts))
        genes_transcripts=map(list,genes_transcripts)
        item.append(genes_transcripts)
        genes_transcripts=[]

    #Appened all transcripts in transcripts list for each gene with the same annotation
    for List in list_of_annotations:
        for item in List[1]:
            for i in object:
                if List[0] == i.Consequence and item[0] == i.Gene_Symbol:
                    transcripts.append(i.Feature_type)

            transcripts=set(transcripts)
            transcripts=list(transcripts)
            item.insert(1, transcripts)
            transcripts = []
        List[0]=List[0].replace('_',' ')
        List[0]=List[0].replace('variant','')



    return list_of_annotations



def Display_variants(ID_list ):

    Variants = []
    filtered_variants = 0
    # Gene_list = []
    GENES_ID_list = []
    for i in ID_list:
        var_list = []
        cal_freq = ''
        # Loop through Unique ID's to Retrieve Variant and Population Records
        variant = get_object_or_404(Variant, id=i)
        population = get_object_or_404(Population, Variant_id=i)

        # Retrieves 1st Annotation of each Annotation list
        Annotations = Annotation.objects.filter(Variant_id=i)

        # Calculates Frequency for each Variant in All Population
        if population.AN_Adjusted != 0:
            cal_freq = float(population.AC_Adjusted) / float(population.AN_Adjusted)
            cal_freq=format(cal_freq,'.4g')

        # Counts Filtered Variants
        if variant.Filter == 'PASS':
            filtered_variants += 1

        # Get the Required Annotations for display in HTML Template

        var_list.append(variant)
        var_list.append(variant.Chromosome_Number)
        var_list.append(variant.Position)
        var_list.append(variant.Filter)
        var_list.append(Annotations[0].HGVS_DNA)
        var_list.append(Annotations[0].Consequence)
        var_list.append(Annotations[0].LoF_flags)
        var_list.append(population.AC_Adjusted)
        var_list.append(population.AN_Adjusted)
        var_list.append(population.AC_Homozygous)
        var_list.append(cal_freq)
        Variants.append(var_list)

        # Retrieves Unique Genes ID for Display in HTML Template
        for ann in Annotations:
            if ann.Gene != "":
                # Gene_list.append(ann.Gene_Symbol)
                GENES_ID_list.append(ann.Gene)

    Genes = set(GENES_ID_list)
    # Genes_ID = set(GENES_ID_list)
    return Variants, filtered_variants, Genes



def calc_allele_freq(object):
    allele_freq=[]
    allele_freq.append(float(object.AC_South_Asian ) / float(object.AN_South_Asian))
    allele_freq.append(float(object.AC_African_American) / float(object.AN_African_American))
    allele_freq.append(float(object.AC_East_Asian) / float(object.AN_East_Asian ))
    allele_freq.append(float(object.AC_Finnish) / float(object.AN_Finnish))
    allele_freq.append(float(object.AC_Non_Finnish_European ) / float(object.AN_Non_Finnish))
    allele_freq.append(float(object.AC_American) / float(object.AN_American ))
    allele_freq.append(float(object.AC_Other) / float(object.AN_Other))
    allele_freq.append(float(object.AC_Adjusted) / float(object.AN_Adjusted))

    # for displaying 4 digits(1-9) only after .
    for i in range(len(allele_freq)):
        allele_freq[i]=format(allele_freq[i],'.4g')

    return allele_freq


    #return HttpResponse("<h1>Details for Gene: " + gene_symbol + "</h1>")

# def ComplexQuery_view(request):
#     if request.is_ajax():
#         if request.method == 'POST':
#             json_string = request.read()
#             dict = json.loads(json_string)
#             print("############################")
#             print(dict)
#             print(json_string)
#             print(dict['query'][1])
#             query_list = Parse(dict)
#             print(query_list)
#             print(query_list[0]['column_name'])
#             return HttpResponse(query_list[1]['operation'])
@csrf_exempt
def ComplexQuery_view(request):
    start_time = time.time();
    ID_list = []
    Transcript_list = []
    if request.is_ajax():
        if request.method == 'POST':
            json_string = request.read()
            dict = json.loads(json_string)
            #json_string = json.loads(request.body)
            #test_var = str(gene_symbol)

            # Check if Request is Gene symbol or Gene_code to Determine Filter Type
            #try:
            #if test_var[:4] == 'ENSG':
            #    Annotations = Annotation.objects.filter(Complex_query(0))
            #    gene_symbol = Annotations[0].Gene_Symbol
            #else:
            #    Annotations = Annotation.objects.filter(Complex_query(0))

            ##test_end
            # print attribute_exists("Chromosome_Number")
            # print attribute_exists("AC_African_American")
            # print attribute_exists("Consequence")
            # print attribute_exists("Mohamad")
            ##test_end

            ##test_start
            # q1 = Q(**{"Variant__" + "Reference": "G"}) & Q(**{"Gene_Symbol": "NRIP1"})  # test
            # Annotations = Annotation.objects.filter(q1)  # test
            ##test_end

            ##q1 = Q(**{"Population__" + "Allele_Count": "29"}) #& Q(**{"Gene_Symbol": "NRIP1"})  # test
            ##Annotations = Annotation.objects.filter(q1)

            #Annotations = Annotation.objects.filter(Complex_query(0, Parse(Query)))

            for Query in Q_list:

                Qobj = Complex_query(0, Parse(Query));

                Annotations = Annotation.objects.filter(Qobj)

                print Annotations.query

            # Annotations = Annotation.objects.filter(Gene_Symbol=gene_symbol)
            gene_ID = Annotations[0].Gene
            gene_symbol=Annotations[0].Gene_Symbol
            # Retrieves Unique Foreign Keys & Transcripts of each GENE
            for rec in Annotations:
                ID_list.append(rec.Variant_id)
                Transcript_list.append(rec.Feature_type)

            Transcripts = set(Transcript_list)
            var_id_list = set(ID_list)


            print '\n'.join(str(p) for p in var_id_list)

            Variants, filtered_variants, Genes_ID = Display_variants(var_id_list)
            query_time = format(float(time.time() - start_time), '.4g')

            # Saving time of query to a file
            #target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
            #target.write("ComplexQuery=")
            #target.write(test_var)
            #target.write("  ")
            #target.write(query_time)
            #target.write("\n")
            #target.close()

            return render(request, 'ComplexQuery_view.html',
                          {'gene_ID': gene_ID, 'Variants_num': len(var_id_list), 'gene_symbol': gene_symbol,
                           'Variants': Variants, 'filtered_variants': filtered_variants, 'Transcripts': Transcripts,
                           'query_time': query_time})
            #except:
            #    return error_view(request)
    else:
        return request


@csrf_exempt
def ComplexQuerydebug_view(request):
    start_time = time.time();
    ID_list = []
    Transcript_list = []

    if 'Query' in request.POST:
        json_string = request.POST['Query'].split('=', 1)[-1]

        #json_string = request.read()
        dict = json.loads(json_string)
        #json_string = json.loads(request.body)
        #test_var = str(gene_symbol)

        # Check if Request is Gene symbol or Gene_code to Determine Filter Type
        #try:
        #if test_var[:4] == 'ENSG':
        #    Annotations = Annotation.objects.filter(Complex_query(0))
        #    gene_symbol = Annotations[0].Gene_Symbol
        #else:
        #    Annotations = Annotation.objects.filter(Complex_query(0))

        ##test_end
        # print attribute_exists("Chromosome_Number")
        # print attribute_exists("AC_African_American")
        # print attribute_exists("Consequence")
        # print attribute_exists("Mohamad")
        ##test_end

        ##test_start
        # q1 = Q(**{"Variant__" + "Reference": "G"}) & Q(**{"Gene_Symbol": "NRIP1"})  # test
        # Annotations = Annotation.objects.filter(q1)  # test
        ##test_end

        ##q1 = Q(**{"Population__" + "Allele_Count": "29"}) #& Q(**{"Gene_Symbol": "NRIP1"})  # test
        ##Annotations = Annotation.objects.filter(q1)

        #Annotations = Annotation.objects.filter(Complex_query(0, Parse(Query)))

        # for Query in Q_list:
        #
        #     Qobj = Complex_query(0, Parse(Query));
        #
        #     Annotations = Annotation.objects.filter(Qobj)
        #
        #     print Annotations.query

        Annotations = Annotation.objects.filter(Complex_query(0, Parse(dict)))

        # Annotations = Annotation.objects.filter(Gene_Symbol=gene_symbol)
        gene_ID = Annotations[0].Gene
        gene_symbol=Annotations[0].Gene_Symbol
        # Retrieves Unique Foreign Keys & Transcripts of each GENE
        for rec in Annotations:
            ID_list.append(rec.Variant_id)
            Transcript_list.append(rec.Feature_type)

        Transcripts = set(Transcript_list)
        var_id_list = set(ID_list)

        print len(var_id_list)

        var_num = len(var_id_list)

        # start_print_time = time.time();
        # print '\n'.join(str(p) for p in var_id_list)
        # print_time = format(float(time.time() - start_print_time), '.4g')

        # print print_time

        Variants, filtered_variants, Genes_ID = Display_variants(var_id_list)
        query_time = format(float(time.time() - start_time), '.4g')

        # Saving time of query to a file
        #target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
        #target.write("ComplexQuery=")
        #target.write(test_var)
        #target.write("  ")
        #target.write(query_time)
        #target.write("\n")
        #target.close()

        return render(request, 'ComplexQuerydebug_view.html',
                      {'gene_ID': gene_ID, 'Variants_num': len(var_id_list), 'gene_symbol': gene_symbol,
                       'Variants': Variants, 'filtered_variants': filtered_variants, 'Transcripts': Transcripts,
                       'query_time': query_time, 'var_num': var_num})
        #except:
        #    return error_view(request)