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