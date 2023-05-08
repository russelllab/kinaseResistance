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
              x: ['Activating', 'Deactivating', 'Resistance', 'Neutral'],
              y: [response['activating'],
                  response['deactivating'],
                  response['resistant'],
                  // response['neutral']
                ],
              marker:{
                color: ['rgba(0, 158, 115, 1.0)',
                        'rgba(213, 94, 0, 1.0)',
                        'rgba(0, 114, 178, 1.0)',
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
              title: 'Predicted effect of '+mutation+'<br>on '+kinase+' activity',
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