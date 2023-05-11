
function makeDom(rows, head, filePath, firstGroup, secondGroup){
  function unpack(rows, key, mutType) {
    return rows.map(function(row) {
                  // alert(row['mutType']+row[key]);
                  if (row['mutType'] == mutType) {
                                return row[key]; 
                  }
                  else {
                                return null;
                  }
                });
  }
  var data = [{
  type: 'violin',
  x: unpack(rows, head, firstGroup),
  y: unpack(rows, 'score', firstGroup),
  legendgroup: firstGroup,
  scalegroup: firstGroup,
  name: firstGroup,
  side: 'negative',
  box: {
    visible: true
  },
  line: {
    color: 'green',
    width: 2
  },
  meanline: {
    visible: true
  }
}, {
  type: 'violin',
  x: unpack(rows, head, secondGroup),
  y: unpack(rows, 'score', secondGroup),
  legendgroup: secondGroup,
  scalegroup: secondGroup,
  name: secondGroup,
  side: 'positive',
  box: {
    visible: true
  },
  line: {
    color: 'red',
    width: 2
  },
  meanline: {
    visible: true
  }
}]

var layout = {
  title: firstGroup + " and " + secondGroup,
  yaxis: {
    zeroline: false
  },
  violingap: 0,
  violingroupgap: 0,
  violinmode: "overlay",
}

return [data, layout];

}

function showFeatures(objects)
{
  d3.tsv(objects['filePath'], function(err, rows){
    [data, layout] = makeDom(rows, objects['head'], objects['filePath'], 
                            objects['firstGroup'], objects['secondGroup']);
    Plotly.newPlot(objects['domElement'] , data, layout);
  });
}
