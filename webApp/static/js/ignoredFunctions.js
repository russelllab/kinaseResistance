/* Define the ignored datatable */
function defineIgnoredDataTable (tableID, uniqID)
{
    // alert(uniqID+tableID);
    $(document).ready(function () {
        var table = $('#'+tableID).DataTable({
            ajax: '/static/predictor/output/'+uniqID+'/ignored.json',
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
                    { data: 'name' }, // input
                    { data: 'reason' }
                ]
            });
    });
}