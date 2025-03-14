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
                            // VerdictNDA
                            if ( data.verdictNDA.includes('Activating') ) {
                                if ( data.verdictNDA.includes('Low') ) {
                                    $('td:eq(5)', row).css('background-color', '#90EE90');
                                }
                                else if ( data.verdictNDA.includes('Medium') ) {
                                    $('td:eq(5)', row).css('background-color', '#32CD32');
                                }
                                else {
                                    $('td:eq(5)', row).css('background-color', '#008000');
                                    $('td:eq(5)', row).css('color', 'white');
                                }

                            }
                            else if  ( data.verdictNDA.includes('Deactivating') ) {
                                if ( data.verdictNDA.includes('Low') ) {
                                    $('td:eq(5)', row).css('background-color', '#FF8A8A');
                                }
                                else if ( data.verdictNDA.includes('Medium') ) {
                                    $('td:eq(5)', row).css('background-color', '#FF5C5C');
                                }
                                else {
                                    $('td:eq(5)', row).css('background-color', '#960018');
                                    $('td:eq(5)', row).css('color', 'white');
                                }
                            }

                            // VerdictR
                            if ( data.verdictR.includes('Resistance') ) {
                                $('td:eq(6)', row).css('background-color', '#0072b2');
                                $('td:eq(6)', row).css('color', 'white');
                            }

                            // ActvDeact
                            if ( data.AIvLD >= 0.5 ) {
                                $('td:eq(8)', row).css('background-color', '#009e73');
                            }
                            else if  ( data.AIvLD < 0.5 ) {
                                $('td:eq(8)', row).css('background-color', '#d55e00');
                            }
                            
                            // NeutralvsDeactvsAct
                            if (data.N == 'NA') {}
                            else if ( data.N >= data.D && data.N >= data.A ) {
                                $('td:eq(11)', row).css('background-color', '#F2E34C');
                            }
                            else if ( data.D >= data.N && data.D >= data.A ) {
                                $('td:eq(10)', row).css('background-color', '#d55e00');
                            }
                            else if ( data.A >= data.D && data.A >= data.N ) {
                                $('td:eq(9)', row).css('background-color', '#009e73');
                            }

                            // ResistantvsNeutral
                            if ( data.RvN >= 0.5 ) {
                                $('td:eq(12)', row).css('background-color', '#0072b2');
                            }
                            
                        },
            dom: 'Bfrtip',
            // buttons: [
            //     'copy', 'csv', 'excel', 'pdf', 'print', '<h3>columnsToggle</h3>',
            //     '<h3 data-title="Information" data-intro="Hover on these buttons to know more about the column.">colvis</h3>'
            // ],
            buttons: [
                {
                    extend: 'collection',
                    text: '<span data-title="Export table" data-intro="Select one of the options here to export the table.">Export</span>',
                    className: 'custom-html-collection',
                    buttons: [
                        'copy',
                        'pdf',
                        'csv',
                        'excel',
                        'print',
                    ]
                },
                {
                    extend: 'collection',
                    text: '<span data-title="Colvis" data-intro="Show/hide columns in the table. By default, only relevant columns are displayed.">Column Visibility</span>',
                    className: 'custom-html-collection',
                    buttons: [
                        'columnsToggle'
                    ]
                }
            ],
            columns: [
                    // {
                    // className: 'dt-control',
                    // orderable: false,
                    // data: null,
                    // defaultContent: '',
                    // },
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
                    { data: 'verdictNDA' },
                    { data: 'verdictR' },
                    { data: 'AIvLD' },
                    { data: 'A' },
                    { data: 'D' },
                    { data: 'N' },
                    // { data: 'AIvNLD' },
                    // { data: 'LDvNAI' },
                    // { data: 'AIvN' },
                    // { data: 'LDvN' },
                    { data: 'RvN' }
                    
                ],
            columnDefs: [
                {
                    targets: [1,3,5,7,11,12,13,14,15],
                    visible: false,
                    searchable: true
                }
            ],
            });

        // Activate this when you want to have a child row
        /*
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
        */
    });
}

function overEntriesAlert(overEntries)
{
    // get a reference to the div element
		myDiv = document.getElementById('overEntriesAlert');
		// check conditions and activate the div if they are met
        // alert (error);

		if (overEntries != 0) {
		myDiv.style.display = 'block';
		} else {
		myDiv.style.display = 'none';
		}
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
    var table_contents = [
                        // 'Info',
                        'View details',
                        'Input',
                        'Gene',
                        'UniProt<br>accession',
                        'Variant',
                        "Site<br>(+/-5 residues)",
                        // 'HMM<br>position',
                        // 'Alignment<br>position',
                        'Kinase domain<br>region',
                        'PTM<br>type',
                        'Known functional<br>consequence',
                        'Verdict<br>(Activating/<br>Deactivating)',
                        'Verdict<br>(Resistance)',
                        'Pred<br>(A vs D)',
                        'Pred A',
                        'Pred D',
                        'Pred N',
                        'Pred R'
                        ];
    
    var table_contents_text = [
                                // 'Known information about the input',
                                'View details of the prediction for the input',
                                'Input',
                                'Gene name',
                                'UniProt accession',
                                'Variant',
                                'Adjacent residues to variant site',
                                // 'hidden Markov Model position',
                                // 'Alignment position',
                                'Location of the variant site in the canonical kinase domain',
                                'Known PTMs at the variant site',
                                'Known functional consequence of the variant',
                                'Likelihood of the variant to be activating/deactivating',
                                'Likelihood of the variant to be resistant',
                                'Predicted probability of the variant to be activating(>=0.5) vs deactivating(<0.5)',
                                'Predicted probability of the variant to be activating',
                                'Predicted probability of the variant to be deactivating',
                                'Predicted probability of the variant to be neutral',
                                'Predicted probability of the variant to be resistance',
                                ];
    
    var table_head_foot_id = ['table-head', 'table-foot'];
    for (let i = 0; i <= 1; i++) {
      var tableText = document.getElementById(table_head_foot_id[i]);
      // Create a new row element
      const row = document.createElement("tr");
      if (i == 0) {
        row.setAttribute("data-title", 'Header');
        row.setAttribute("data-intro", 'Description of the table columns.');
        row.setAttribute("style", "text-align: center;");
      }
      
      // Create 3 cells in the row using a nested for loop
      for (let j = 0; j < table_contents.length; j++) {
        // Create a new cell element
        const cell = document.createElement("th");

        // const cellText = document.createTextNode(table_contents[j]);

        // cell.appendChild(cellText);
        // cell.innerHTML = table_contents[j] + '<br>';
        cell.innerHTML = table_contents[j];
        cell.style = "white-space: nowrap";

        const but = document.createElement("span");
        but.setAttribute("data-toggle", "tooltip");
        but.setAttribute("data-placement", "top");
        but.setAttribute("title", table_contents_text[j]);
        but.style = "display: inline-block; padding-left: 5px; padding-right: 5px";

        // // Information icon
        // const icon = document.createElement("i");
        // icon.setAttribute("class", "bi bi-info-circle");
        // if ((i == 0) && (j == 0)) {
        //     icon.setAttribute("data-title", 'Information');
        //     icon.setAttribute("data-intro", 'Hover on these buttons to know more about the column.');
        //   }

        // but.appendChild(icon);
        cell.appendChild(but);
                
        // Add the cell to the row
        row.appendChild(cell);
      }
      
      // Add the row to the table body
      tableText.appendChild(row);
    
    }
}