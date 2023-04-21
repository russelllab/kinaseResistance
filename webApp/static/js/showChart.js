function showChart(uniqID, kinase, mutation, results)
{
    $.ajax({
		url: '/AJAXChart',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({'uniqID': uniqID, 'kinase': kinase, 'mutation': mutation, 'results': results}),
		success: function (response){
            var trace1 = {
              x: ['Activating', 'Deactivating', 'Drug resistant', 'Neutral'],
              y: [response['activating'],
                  response['deactivating'],
                  response['resistant'],
                  response['neutral']
                ],
              marker:{
                color: ['rgba(50,171, 96, 0.7)',
                        'rgba(222,45,38,0.8)',
                        'rgba(204,204,204,1)',
                        'rgba(204,204,204,1)'
                      ]
              },
              type: 'bar'
            };
            
            var data = [trace1];
            var layout = {
              title: 'Predicted effect of the mutation on the kinase activity',
            };

            Plotly.newPlot('myChart', data, layout);

		}
	});
}