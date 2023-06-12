function showChart(uniqID, kinase, mutation, results)
{
    $.ajax({
		url: '/AJAXChart',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({'uniqID': uniqID, 'kinase': kinase, 'mutation': mutation, 'results': results}),
		success: function (response){
            yValues = [response['activating'],
                      response['deactivating'],
                      response['resistance'],
                      // response['activating_AIvLD'],
                      // response['deactivating_AIvLD']
                      ];
            var trace1 = {
              name : 'Prediction',
              x: ['Activating', 'Deactivating', 'Resistance'],
              y: yValues,
              marker:{
                color: ['rgba(0, 158, 115, 1.0)',
                        'rgba(213, 94, 0, 1.0)',
                        'rgba(0, 114, 178, 1.0)',
                        // 'rgba(0, 158, 115, 1.0)',
                        // 'rgba(213, 94, 0, 1.0)',
                      ]
              },
              type: 'bar',
              text: yValues.map(String),
              textposition: 'auto',
            };
            
            var data = [trace1];
            var layout = {
              title: 'Effect of '+mutation+'<br>on '+kinase+' activity',
              yaxis: {
                // autorange: true,
                range: [0, 1.0],
                type: 'linear'            
              }
            };

            Plotly.newPlot('predictionChart1', data, layout);

            yValues = [
                      // response['activating_AIvLD'].toFixed(3),
                      response['activating_AIvLD'],
                      ];
            var trace1 = {
              name : 'Activating',
              x: ['A vs D'],
              y: yValues,
              marker:{
                color: [
                        'rgba(0, 158, 115, 1.0)',
                      ]
              },
              type: 'bar',
              text: yValues.map(String),
              textposition: 'auto',
            };

            yValues = [
              response['deactivating_AIvLD']
              ];
            var trace2 = {
              name : 'Deactivating',
              x: ['A vs D'],
              y: yValues,
              marker:{
                color: [
                        'rgba(213, 94, 0, 1.0)',
                      ]
              },
              type: 'bar',
              text: yValues.map(String),
              textposition: 'auto',
            };
            
            var data = [trace1, trace2];
            var layout = {
              title: 'Effect of '+mutation+'<br>on '+kinase+' activity',
              yaxis: {
                // autorange: true,
                range: [0, 1.0],
                type: 'linear'            
              },
              barmode: 'stack'
            };

            Plotly.newPlot('predictionChart2', data, layout);

		}
	});
}