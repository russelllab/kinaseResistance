function format(d) {
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
    return (
        '<table cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">' +
        '<tr>' +
        d.text +
        // textToDisplay +
        '</tr>' +
        '</table>'
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
                            // alert(data.prediction)
                            // color based on score
                            if ( data.predAD >= 0.5 ) {
                                // $(row).addClass( 'important' );
                                // alert(data.predAD);
                                $('td:eq(9)', row).css('background-color', 'Lightgreen');
                            }
                            else if ( data.predAD < 0.5 ) {
                                $('td:eq(9)', row).css('background-color', '#FF6863');
                            }
                            if ( data.predRN >= 0.5 ) {
                                // $(row).addClass( 'important' );
                                // alert(data.predRN);
                                $('td:eq(10)', row).css('background-color', '#7CB9E8');
                            }
                            else if ( data.predRN < 0.5 ) {
                                // $('td:eq(10)', row).css('background-color', '#FF6863');
                            }
                        },
            dom: 'Bfrtip',
            buttons: [
                'copy', 'csv', 'excel', 'pdf', 'print'
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
                    { data: 'hmmPos' },
                    { data: 'ptmType' },
                    { data: 'mutType' },
                    { data: 'predAD' },
                    { data: 'predRN' }
                ]
            });
    
        $('#'+tableID+' tbody').on('click', 'td.dt-control', function () {
            var tr = $(this).closest('tr');
            var row = table.row(tr);

            if (row.child.isShown()) {
                // This row is already open - close it
                row.child.hide();
                tr.removeClass('shown');
            } else {
                // Open this row
                row.child(format(row.data())).show();
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