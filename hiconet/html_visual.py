HTML_HEAD = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"> 
        <html xmlns="http://www.w3.org/1999/xhtml"> 
        <head> 
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" /> 
        <title>HiCoNet Report</title> 

        <style type="text/css"> 
        body {width: 960px; padding: 10px; }
        
        div.colorbar {margin-bottom: 7px}
        div.networkvisual { margin: auto; height: 600px; width: 960px; overflow: auto !important;}
        div.svg2 { margin: auto; height: 600px; width: 960px; overflow: auto !important;}

        div.stats {margin-top: 10; font-size:0.7em; color:#888;}

        div.moduleline{font-size: 0.9em; font-weight: bold;
        }
        div.metabolites {
            font-size: 0.7em;
            padding-left: 15px; padding-bottom: 7px;
        }
        
        footer { margin-top: 30px;    padding: 10px 0;    border-top: 1px solid #582E2E;    font-size:0.7em;    color:#888;}
        
        h1, h2, h3, h4, h5, h6 {
            font-family: 'Trebuchet MS', 'Lucida Grande', Arial, Sans-Serif;
            font-weight: bold;
        }
        
        h1 { font-size: 1.4em; }
        h2 { font-size: 1.2em; }
        h3 { font-size: 1em; }
        h4 { font-size: 0.9em;     padding-left: 15px; padding-bottom: 2px; margin-bottom:0px;}
        th {
            background: #DAFDCF;
            font-size:0.9em;
            text-align: left;
            padding:5px;
        }
        
        tr:nth-child(even) {background-color: #f2f2f2}
        tr:hover {background-color: #f5f5f5}
        
        td {
            font-family: Verdana, Arial, sans-serif;
            color: black;
            font-size:0.7em;
            margin:10px;
            margin-top:0px;
            padding: 5px;
        }
        
        .node {
          stroke: #fff;
          stroke-width: 1.5px;
        }
        
        .link {
          stroke: #999;
          stroke-opacity: .6;
        }
        div.network_selection {float: left; margin-top: 6; margin-bottom: 10; font-size:0.8em;}
        </style> 
        <script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>
        <script src="http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js" charset="utf-8"></script>
        <script src="http://mummichog.org/download/cytoscape.min.js" charset="utf-8"></script>

        </head> 
        <body> 
        <h1>cytoscape.js play</h1>



        <div id="networkvisual" class="networkvisual" ></div>
        """
        
HTML_END = """</body> </html> """
        
javascript_HEAD = """
        <script type="text/javascript" charset="utf-8">
        """
        
javascript_END = """
        var w = 960, h = 600;
        var color = d3.scale.category20b();
        var node_sizes = [16, 8, 32]
        
        var scalebar = d3.select("#colorbar").append("svg").attr("id", "svg_colorbar").attr("width", w).attr("height", 10);
        var coloridx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19];
        scalebar.selectAll("rect")
           .data(coloridx)
           .enter()
           .append("rect")
           .attr("x", function(d, i) {return i * 48;})
           .attr("y", 0)
           .attr("width", 48)
           .attr("height", 10)
           .attr("fill", function(d, i) { return color(i); });
           
        
        // initial draw
        cyto_draw_figure(cytonodes[0], cytoedges[0], node_sizes[0]);
        var cy = $('#networkvisual').cytoscape('get');
        cy.panBy({ x: 0, y: -200 });
        
        // support cytoscape.js  
        function cyto_draw_figure(nodes, links, node_size) {  

                    var svg = document.createElement('div');
                    svg.setAttribute("id","svg2");
                    svg.setAttribute("class","svg2");
                    var nw = document.getElementById("networkvisual");
                    nw.appendChild(svg);

                $('#svg2').cytoscape({
                        style: cytoscape.stylesheet()
                        .selector('node')
                          .css({
                            'content': 'data(id)',
                            'font-size': node_size,
                            'min-zoomed-font-size': 8,
                            'background-color': 'mapData(group, 1, 2, blue, red)'
                          })
                        .selector('edge')
                          .css({
                            'target-arrow-shape': 'none',
                            'width': 4,
                            'line-color': '#ddd',
                            'target-arrow-color': '#ddd'
                          })
                        .selector('.highlighted')
                          .css({
                            'background-color': '#61bffc',
                            'line-color': '#61bffc',
                            'target-arrow-color': '#61bffc',
                            'transition-property': 'background-color, line-color, target-arrow-color',
                            'transition-duration': '0.5s'
                          }),
  
                      layout: {
                        name: 'concentric',
                        minDist: 40,
                        directed: true,
                        padding: 10,
                        //fit: false,
                      },
                      
                      elements: {
                          nodes: nodes, 
                          edges: links
                        },
  
                      zoom: 1,
                      minZoom: 0.1,
                       maxZoom: 5,

                });
    
        }
    
    
        // added support of two styles
        function updateDrawing() {
            cyto_draw_figure(cytonodes[0], cytoedges[0], node_sizes[0]);
        }
    
    </script>
    """


def quote(s):
    return '"'+str(s)+'"'

def make_cytoscape_js_from_json(json_input):
    pass

def make_js_from_network(combined_network):
    '''
    Make code for cytoscape.js visualization.
    
    Note: combined_network is a list of edges: [(network_name, g, m, PLSscore, p-value ), ...]
    
    cytoscape.js visualization:
        var cytonodes = [ [
            { data: { id: "C01227", weight: 2.11 } },
            { data: { id: "C05475", weight: 3.65 } },
            { data: { id: "C05487", weight: 3.65 } },
            { data: { id: "C03205", weight: -1.8 } },
                  ], ... ]
  
                var cytoedges = [ [
            { data: { id: "C01227C05498", weight: 1, source: "C01227", target: "C05498" } },
            { data: { id: "C01227C05487", weight: 1, source: "C01227", target: "C05487" } },
            { data: { id: "C01227C03205", weight: 1, source: "C01227", target: "C03205" } }, ... ]]
    
    Color coded by group. Thru additional argument colordict.

    '''
    cynodestr, cyedgestr = '[', '['
    nodes_society_1, nodes_society_2 = set([e[1] for e in combined_network]), set([e[2] for e in combined_network])
    for n in nodes_society_1:
        cynodestr += '{ data: ' + '{id:%s, group:1} }, ' %quote(n)
    for n in nodes_society_2:
        cynodestr += '{ data: ' + '{id:%s, group:2} }, ' %quote(n)
    for e in combined_network:
        # [(network_name, g, m, PLSscore, p-value ), ...]
        e = [str(x) for x in e]
        # tweaking str types
        cyedgestr += '{ data: { id: "%s", weight: %s, source: %s, target: %s } }, ' %('-'.join(e[1:3]), e[3], quote(e[1]), quote(e[2]))

    cynodestr += '], '
    cyedgestr += '], '
    
    total_cytoscapejs_data = '        var cytonodes = [ ' + cynodestr + '];\n\n        var cytoedges = [' + cyedgestr + '];\n\n'
    return total_cytoscapejs_data
    

def make_html_page(outfile, js_data):
    # gloabl HTML_HEAD, HTML_END, javascript_HEAD, javascript_END
    with open(outfile, "w") as O:
        O.write(HTML_HEAD + javascript_HEAD + js_data + javascript_END + HTML_END)


