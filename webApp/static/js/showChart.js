function showChart(uniqID, kinase, mutation, results)
{
    $.ajax({
		url: '/AJAXChart',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({'uniqID': uniqID, 'kinase': kinase, 'mutation': mutation, 'results': results}),
		success: function (response){
            // console.log(response);
            // AIvLD preditcion
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
              // title: 'Effect of '+mutation+'<br>on '+kinase+' activity',
              title: 'Prediction of A vs D',
              yaxis: {
                // autorange: true,
                range: [0, 1.0],
                type: 'linear'            
              },
              barmode: 'stack'
            };

            Plotly.newPlot('predictionChart1', data, layout);
            // NvDvA preditcion
            yValues = [
                      response['activating'],
                      response['deactivating'],
                      response['neutral'],
                      // response['activating_AIvLD'],
                      // response['deactivating_AIvLD']
                      ];
            var trace1 = {
              name : 'Prediction',
              x: ['Activating', 'Deactivating', 'Neutral'],
              y: yValues,
              marker:{
                color: [
                        'rgba(0, 158, 115, 1.0)',
                        'rgba(213, 94, 0, 1.0)',
                        'rgba(242, 227, 76, 1.0)',
                        
                        // 'rgba(0, 114, 178, 1.0)',
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
              // title: 'Effect of '+mutation+'<br>on '+kinase+' activity',
              title: 'Prediction of N vs D vs A',
              yaxis: {
                // autorange: true,
                range: [0, 1.0],
                type: 'linear'            
              }
            };

            Plotly.newPlot('predictionChart2', data, layout);

            // RvN preditcions
            yValues = [
              // response['activating_AIvLD'].toFixed(3),
              response['resistance'],
              ];
            var trace1 = {
              name : 'Resistance',
              x: ['R vs N'],
              y: yValues,
              marker:{
                color: [
                        'rgba(0, 114, 178, 1.0)',
                      ]
              },
              type: 'bar',
              text: yValues.map(String),
              textposition: 'auto',
            };

            // yValues = [
            //   1 - response['resistance'],
            //   ];
            // var trace2 = {
            //   name : 'Neutral',
            //   x: ['R vs N'],
            //   y: yValues,
            //   marker:{
            //     color: [
            //             'rgba(242, 227, 76, 1.0)',
            //           ]
            //   },
            //   type: 'bar',
            //   text: yValues.map(String),
            //   textposition: 'auto',
            // };
            
            // var data = [trace1, trace2];
            var data = [trace1];
            var layout = {
              // title: 'Effect of '+mutation+'<br>on '+kinase+' activity',
              // title: 'Prediction of R vs N',
              title: 'Prediction of R',
              yaxis: {
                // autorange: true,
                range: [0, 1.0],
                type: 'linear'            
              },
              barmode: 'stack'
            };

            Plotly.newPlot('predictionChart3', data, layout);

		}
	});
}