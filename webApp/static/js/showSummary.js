function showSummary(uniqID, kinase, mutation, results)
{
    $.ajax({
		url: '/AJAXSummary',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({'uniqID': uniqID, 'kinase': kinase, 'mutation': mutation, 'results': results}),
		success: function (response){
            document.getElementById('summaryCard').innerHTML = response['text'];

		}
	});
}

function defineSummaryDataTable (tableID, uniqID, kinase, mutation)
{	
    $.ajax({
		url: '/AJAXSummaryTable',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({'uniqID': uniqID, 'kinase': kinase, 'mutation': mutation}),
		success: function (response){
					$(document).ready(function () {
						var table = $('#'+tableID).DataTable({
							// scrollY:        "300px",
							// scrollX:        true,
							// scrollCollapse: true,
							// paging:         false,

							data: response['data'],
							dom: 'Bfrtip',
							lengthChange: false,
							buttons: [
								'copy', 'csv', 'excel', 'pdf', 'print', 'colvis',
								{
									text: "Show all instances",
									action: function(e, dt, node, config){
										dt.column(7).search('').draw();
										dt.column(1).search('').draw();
									}
								},
								{
									text: "Show only "+response['gene']+" instances",
									action: function(e, dt, node, config){
										dt.column(7).search('').draw();
										dt.column(1).search(response['acc']).draw();
									}
								},
								{
									text: "Show only activating instances",
									action: function(e, dt, node, config){
										dt.column(1).search('').draw();
										dt.column(7).search("^"+'activating', true, false).draw();
									}
								},
								{
									text: "Show only resistance instances",
									action: function(e, dt, node, config){
										dt.column(1).search('').draw();
										dt.column(7).search("^"+'resistance').draw();
									}
								},
								{
									text: "Show only deactivating instances",
									action: function(e, dt, node, config){
										dt.column(1).search('').draw();
										dt.column(7).search("^"+'deactivating').draw();
									}
								},
							],
							columnDefs: [
								{ width: 200, targets: 7 },
								{
									targets: "_all",
									className: 'dt-center'
								}
							],
							columns: [
								{ title: 'Gene<br>Name' },
								{ title: 'UniProt<br>Acc' },
								{ title: 'WT' },
								{ title: 'Position' },
								{ title: 'MUT' },
								{ title: 'HMM<br>position' },
								{ title: 'Alignment<br>position' },
								{ title: 'Info<br>type' },
								{ title: 'Description' },
								{ title: 'Reference' },
								],
							});
						table.buttons().container()
						.appendTo( '#example_wrapper .col-md-6:eq(0)' );
					});
				}
	});
}

function makeSummaryTableHeadFoot()
{   
    var table_contents = ['Info',
                        'Gene<br>name',
                        'UniProt<br>accession',
                        'Sequence<br>position',
                        'HMM<br>position',
                        'Info<br>type',
						'Description',
						'Links'
                        ];
    
    var table_contents_text = ['Known information about the input',
                                'View individual results',
                                'Input mutation',
                                'Gene name',
                                'UniProt accession',
                                'Mutation',
                                'HMM position',
                                'Region of the mutation site<br>in the kinase canonical structure',
                                'Known PTM at the position of the mutation',
                                'Known Activating/Deactivating/Resistance mutation',
                                'Predicted probability of Activating(>=0.5)<br>or Deactivating(<0.5)',
                                'Predicted probability of Resistance'
                                ];
    
    var table_head_foot_id = ['summary-table-head', 'summary-table-foot'];
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