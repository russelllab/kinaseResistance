function showChart(uniqID, kinase, mutation, results)
{
    $.ajax({
		url: '/AJAXChart',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({'uniqID': uniqID, 'kinase': kinase, 'mutation': mutation, 'results': results}),
		success: function (response){
            var trace1 = {
              name : 'Prediction',
              x: ['Activating', 'Deactivating', 'Drug resistant', 'Neutral'],
              y: [response['activating'],
                  response['deactivating'],
                  response['resistant'],
                  // response['neutral']
                ],
              marker:{
                color: ['rgba(50,171, 96, 0.7)',
                        'rgba(222,45,38,0.8)',
                        'rgba(124, 185, 232, 1.0)',
                        // 'rgba(204,204,204,1)'
                      ]
              },
              type: 'bar'
            };
            
            var trace2 = {
              // x: ['Activating', 'Deactivating', 'Drug resistant', 'Neutral'],
              // y: [0.5, 0.5, 0.5],
              // type: 'line'
            };
            var data = [trace1, trace2];
            var layout = {
              title: 'Predicted effect of the mutation on the kinase activity',
              yaxis: {
                // autorange: true,
                range: [0, 1.0],
                type: 'linear'            
              }
            };

            Plotly.newPlot('predictionChart', data, layout);

		}
	});
}