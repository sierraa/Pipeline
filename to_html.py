# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 12:21:09 2015

@author: Sierra Anderson

Generate HTML file displaying comparative analysis results. 

"""

import os
import webbrowser
import string
import generate_text_descriptions

def generate_params(ca_params):
    """ Generate parameters for building html page. 
    ca_params -- parameters for comparative analysis 
    
    """
    
    page_params = dict()
    out = "" 
    if ca_params["output_dir"] != "current":
        out = ca_params["output_dir"] + "/"
    
    if ca_params["enrichment"][0] == "t":
        if ca_params["multiple_comparisons"][0] == "f":
            page_params["enrichment"] = out + ca_params["enrichment_test"] + "_pvals.tab"
        else:
            # need to get file names for corrected pairwise
            page_params["enrichment"] = out + "tukey.tab"
    
    if ca_params["area_plot"][0] == "t" and os.path.isfile(out + "area_plot.png"):
        # extra check for area_plot.png since if there are too many data points
        # the area plot will not be generated 
        page_params["area_plot"] = out + "area_plot.png"
    else:
        ca_params["area_plot"] = "false"
    
    if ca_params["pcoa"][0] == "t":
        page_params["pcoa"] = out + "pcoa_" + ca_params["dist_type"] + ".png"
    
    if ca_params["pca"][0] == "t":
        page_params["pca"] = out + "pca.png"
    
    return page_params
    
def perform_write(ca_params, name):
    """ Helper function to determine if component should be written to html page.
    
    ca_params -- parameters from comparative analysis 
    name -- name of the component to test if it should be written
    
    """
    return ca_params[name][0] == 't'

def write_page(page_params, ca_params, filename='results.html'):   
    """Write an HTML page containing results from comparative analysis script. 
    Open page if specified in ca_params.
    
    page_params -- parameters for building the page 
    ca_params -- parameters from comparative analysis script
    filename -- name of the html file (default='results.html') 
    
    """
    
    contents = ''' <html><head>
        <style>
        body {
            background-color: #D4DCE1;
            font-family: Verdana, sans-serif;
        }
        
        h2 {
            text-align: center;
        }
        
        table.center {
            border: 5px solid #ADBBC3;    
        }
        
        td {
            font-size: 95%;
            padding: 5px;
        }
        
        th {
            font-weight: bold; 
        }
        
        div.img {
            display: block;
            margin-left: auto;
            margin-right: auto;
        }
        
        div.enrich {
            overflow: auto; 
            max-height: 400px;        
        }
        
        img.center {
            display: block;
            margin-left: auto;
            margin-right: auto;
            border: 5px solid #ADBBC3;
        }
        
        table.center {
            margin-left:auto; 
            margin-right:auto;
            background-color: white;
        }
        
        a:link {
            padding: 10px;
            font-weight: bold;
            font-color: black;
            text-decoration: none;            
            background-color: #ADBBC3;
        }
        
        a:visited {
            color: black;
            font-weight: bold;
            text-decoration: none;
        }
        
        a:hover {
            background-color: #D4DCE1;
        }
        
        table.links {           
            border-spacing: 15px;
        }
    
        p.about {
            margin-left: 2cm;
        }        
        
        </style>
        <title>Results</title>
        </head><body>
        <h1>Comparative Analysis Results</h1>
        <p><table style="width:"700" class="links"> '''    
    
    if perform_write(ca_params, 'area_plot'):
        contents += '<th class="links"><a href="#area">Area Plot</a></th>'
            
    if perform_write(ca_params, 'enrichment'):
        contents += '<th class="links"><a href="#enrichment">Enrichment</a></th>'
        
    if perform_write(ca_params, 'pcoa'):
        contents += '<th class="links"><a href="#pcoa">Principal Coordinate Analysis</a></th>'
        
    if perform_write(ca_params, 'pca'):
        contents += '<th class="links"><a href="#pca">Principal Component Analysis</a></th>'
        
    contents += '</table></p>'
        
    contents += '<p>Name: ' + string.capwords(ca_params['name']) + '<br>'
    contents += 'Year: ' + ca_params['year'] + '<br>'
    contents += 'Sequence Type: ' + string.capwords(ca_params['sequence_type']) + '<br>'
    contents += 'Collaborator: ' + string.capwords(ca_params['collaborator']) + '</p>'
    
    page_params = generate_text_descriptions.generate(page_params, ca_params)    
    
    # Area Plot
    if perform_write(ca_params, 'area_plot'):
        contents += '<h2><a name="area">Area Plot</a></h2>'
        contents += '''<p><div class="img">
        <img src="''' + page_params['area_plot'] + '''" width="900" class="center">
        </div></p>'''
        contents += '<p class="about">' + page_params['area_plot_text'] + '</p>'
    
    # PCA
    if perform_write(ca_params, 'pca'):
        contents += '<h2><a name="pca">Principal Component Analysis</a></h2>'
        contents += '''<p><div class="img">
        <img src="''' + page_params['pca'] + '''"width="900" class="center">
        </div></p>'''
    
    # PCoA    
    if perform_write(ca_params, 'pcoa'):
        contents += '<h2><a name="pcoa">Principal Coordinate Analysis</a></h2>'
        contents += '''<p><div class="img">
        <img src="''' + page_params['pcoa'] + '''"width="900" class="center">
        </div></p>'''
        contents += '<p class="about">' + page_params['pcoa_text'] + '</p>'
    
    # Enrichment 
    if perform_write(ca_params, 'enrichment'):
        contents += '<h2><a name="enrichment">Enrichment</a></h2>'
        contents += '<p class="about">' + page_params["enrichment_text"] + '</p>'
        contents += '<div id="enrich"><table style="width:"900" class="center">'
        with open(page_params['enrichment'], 'r') as f:
            i = 0
            for l in f:
                
                if i == 0:
                    l = l.split('\t')
                    contents += '<tr>'
                    for word in l:
                        contents += '<th>' + word.lower().capitalize() + '</th>'
                    contents += '</tr>'
                else:
                    l = l.split('\t')
                    contents += '<tr>'
                    for word in l:
                        contents += '<td>' + word.capitalize() + '</td>'
                    contents += '</tr>'
                    
                i += 1 
            
        contents += '</table></div>'
                
    contents += '</body></html>'
    out = open(filename,'w')    
    out.write(contents)
    out.close()
    
    if ca_params["open_page"][0] == 't':
        webbrowser.open(filename)
