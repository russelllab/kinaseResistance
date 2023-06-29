function format(d, uniqID) {
    // `d` is the original data object for the row
    // var textToDisplay = '';
    // if (d.ptmType != '-') {
    //     textToDisplay += d.name.slice(0,-1) + ' is a known ' + d.ptmType + ' site. <br>'
    // }
    // if (d.mutType != '-') {
    //     var mutationTypes = d.mutType.split('');
    //     for (let i = 0; i<mutationTypes.length; i++) {
    //         if (mutationTypes[i] == 'A') {
    //             mutationTypes[i] = 'activating';
    //         }
    //         else if (mutationTypes[i] == 'D') {
    //             mutationTypes[i] = 'deactivating';
    //         }
    //         else if (mutationTypes[i] == 'R') {
    //             mutationTypes[i] = 'drug-resistant';
    //         }
    //     }
    //     textToDisplay += d.name + ' is a known ' + mutationTypes.join(' and ') + ' mutation. <br>'
    // }
    const name = d.name;
    const kinase = name.split('/')[0];
    const mutation = name.split('/')[1];
    const divElement = document.createElement("div");
    const table = document.createElement("table");
    table.setAttribute("id", "summaryTableOutput");
    table.className = "display order-column";
    table.style = "text-align: center;";
    const header = document.createElement("h5");
    header.innerHTML = "Information known in other kinases";
    header.style = "text-align: left;";
    divElement.appendChild(header);
    divElement.appendChild(table);
    defineSummaryDataTable ("summaryTableOutput", uniqID, kinase, mutation);

    return (
        // '<table id="summaryTable" cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">' +
        // '<table id="summaryTableOutputput" class="display order-column" style="text-align: center;">' +
        // '<tr>' +
        // d.text +
        // textToDisplay +
        // '</tr>' +
        // '</table>'+
        // '<script>'+
        // 'defineSummaryDataTable ("summaryTableOutput", '+uniqID+', '+d.gene+', '+d.mutation+', '+d.acc+');'+
        // '</script>'
        divElement
    );
}


/* Define the datatable */
function defineDataTable (tableID, uniqID)
{
    // alert({{uniqID|safe}});
    $(document).ready(function () {
        var table = $('#'+tableID).DataTable({
            ajax: '/static/predictor/output/'+uniqID+'/output.json',
            // ajax: '/static/predictor/output/GPBXK/output.json',
            "createdRow": function( row, data, dataIndex ) {
                            if (data.N == 'NA') {}
                            
                            else if ( data.N >= data.D && data.N >= data.A ) {
                                $('td:eq(10)', row).css('background-color', '#F2E34C');
                            }
                            else if ( data.D >= data.N && data.D >= data.A ) {
                                $('td:eq(11)', row).css('background-color', '#d55e00');
                            }
                            else if ( data.A >= data.D && data.A >= data.N ) {
                                $('td:eq(12)', row).css('background-color', '#009e73');
                            }
                            if ( data.RvN >= 0.5 ) {
                                $('td:eq(14)', row).css('background-color', '#0072b2');
                            }
                            if ( data.AIvLD >= 0.5 ) {
                                $('td:eq(13)', row).css('background-color', '#009e73');
                            }
                            else if  ( data.AIvLD < 0.5 ) {
                                $('td:eq(13)', row).css('background-color', '#d55e00');
                            }
                        },
            dom: 'Bfrtip',
            buttons: [
                'copy', 'csv', 'excel', 'pdf', 'print', 'colvis'
            ],
            columns: [
                    {
                    className: 'dt-control',
                    orderable: false,
                    data: null,
                    defaultContent: '',
                    },
                    { data: 'view' },
                    { data: 'name' }, // input
                    { data: 'gene' },
                    { data: 'acc' },
                    { data: 'mutation' },
                    { data: 'adjacentSites' },
                    // { data: 'hmmPos' },
                    // { data: 'alnPos' },
                    { data: 'region' },
                    { data: 'ptmType' },
                    { data: 'mutType' },
                    { data: 'N' },
                    { data: 'D' },
                    { data: 'A' },
                    // { data: 'AIvNLD' },
                    // { data: 'LDvNAI' },
                    // { data: 'AIvN' },
                    // { data: 'LDvN' },
                    { data: 'AIvLD' },
                    { data: 'RvN' }
                    
                ]
            });
    
        $('#'+tableID+' tbody').on('click', 'td.dt-control', function () {
            var tr = $(this).closest('tr');
            var row = table.row(tr);
            // alert(row.data().name);

            if (row.child.isShown()) {
                // This row is already open - close it
                row.child.hide();
                tr.removeClass('shown');
            } else {
                // Open this row
                row.child(format(row.data(), uniqID)).show();
                tr.addClass('shown');
            }
        });
    });
}

function ignoredAlert(error)
{
    // get a reference to the div element
		myDiv = document.getElementById('ignoredAlert');
		// check conditions and activate the div if they are met
        // alert (error);
		if (error != "0") {
		myDiv.style.display = 'block';
		} else {
		myDiv.style.display = 'none';
		}

        myDiv = document.getElementById('errorText');
        if (error != "1") {
            myDiv.innerHTML = error + ' inputs were ignored.';
            } else {
                myDiv.innerHTML = error + ' input was ignored.';
            }
}

function makeIgnoredLink(uniqID)
{
    const ignoredLink = document.getElementById('ignored-link');
    // String concatenation
    ignoredLink.href = '/ignored/'+uniqID;

}

function makeTableHeadFoot()
{   
    var table_contents = ['Info',
                        'View',
                        'Input',
                        'Gene<br>name',
                        'UniProt<br>accession',
                        'Mutation',
                        "Site<br>(+/-5 residues)",
                        // 'HMM<br>position',
                        // 'Alignment<br>position',
                        'Region',
                        'PTM<br>type',
                        'Known<br>ADR',
                        'Pred N',
                        'Pred D',
                        'Pred A',
                        'Pred<br>(A vs D)',
                        'Pred R'
                        ];
    
    var table_contents_text = ['Known information about the input',
                                'View individual results',
                                'Input mutation',
                                'Gene name',
                                'UniProt accession',
                                'Mutation',
                                'Adjacent residues to the mutation site',
                                // 'hidden Markov Model position',
                                // 'Alignment position',
                                'Region of the mutation site<br>in the kinase canonical structure',
                                'Known PTM at the position of the mutation',
                                'Known Activating/Deactivating/Resistance mutation',
                                'Predicted probability of Neutral',
                                'Predicted probability of Deactivating',
                                'Predicted probability of Activating',
                                'Predicted probability of Activating(>=0.5) vs Deactivating(<0.5)',
                                'Predicted probability of Resistance',
                                ];
    
    var table_head_foot_id = ['table-head', 'table-foot'];
    for (let i = 0; i <= 1; i++) {
      var tableText = document.getElementById(table_head_foot_id[i]);
      // Create a new row element
      const row = document.createElement("tr");
      if (i == 0) {
        row.setAttribute("data-title", 'Header');
        row.setAttribute("data-intro", 'Description of the table columns.');
      }
      
      // Create 3 cells in the row using a nested for loop
      for (let j = 0; j < table_contents.length; j++) {
        // Create a new cell element
        const cell = document.createElement("th");

        // const cellText = document.createTextNode(table_contents[j]);

        // cell.appendChild(cellText);
        cell.innerHTML = table_contents[j] + '<br>';
        cell.style = "white-space: nowrap";

        const but = document.createElement("span");
        but.setAttribute("data-toggle", "tooltip");
        but.setAttribute("data-placement", "top");
        but.setAttribute("title", table_contents_text[j]);
        but.style = "display: inline-block; padding-left: 5px; padding-right: 5px";

        const icon = document.createElement("i");
        icon.setAttribute("class", "bi bi-info-circle");
        if ((i == 0) && (j == 0)) {
            icon.setAttribute("data-title", 'Help');
            icon.setAttribute("data-intro", 'Click the help buttons to know more.');
          }

        but.appendChild(icon);
        cell.appendChild(but);
                
        // Add the cell to the row
        row.appendChild(cell);
      }
      
      // Add the row to the table body
      tableText.appendChild(row);
    
    }
}