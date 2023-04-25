function showAlignment(uniqID, kinase, mutation, results)
{
    $.ajax({
		url: '/AJAXAlignment',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({'uniqID': uniqID, 'kinase': kinase, 'mutation': mutation, 'results': results}),
		success: function (response){
            // document.getElementById('alignmentCard').innerHTML = response['filepath'];
			const image = document.getElementById('aliView');
      		image.setAttribute('src', response['filepath']);

		}
	});
}