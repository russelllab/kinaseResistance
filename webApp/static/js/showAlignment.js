function openImage(uniqID) {
	// Get the image source
	// const imgSrc = document.getElementById("image-container").src;
	const svgImage = document.getElementById("aliView");
	// Get the SVG image URL
	// const svgImageURL = window.URL.createObjectURL(new Blob([svgImage.outerHTML], {type: "image/svg+xml"}));
	const svgImageURL = 'static/predictor/output/'+uniqID+'/'+svgImage.name;
	alert(svgImageURL);

	// Open the image in a new tab
	window.open(svgImageURL, "_blank");
	}

function downloadSVG(uniqID) {
	// Get the SVG image element
	const svgImage = document.getElementById("aliView");
	// Get the SVG image URL
	// const svgImageURL = window.URL.createObjectURL(new Blob([svgImage.outerHTML], {type: "image/svg+xml"}));
	const svgImageURL = 'static/predictor/output/'+uniqID+'/'+svgImage.name;
	// alert(svgImageURL);
	// Create a link element to download the SVG image
	const link = document.createElement("a");
	link.href = svgImageURL;
	link.download = "image.svg";
	document.body.appendChild(link);
	link.click();
	document.body.removeChild(link);
}
function showAlignment(uniqID, kinase, mutation, results, currentValueWS, currentValueTopN)
{
	if (currentValueWS === undefined) {
		currentValueWS = 30;
	  }
	if (currentValueTopN === undefined) {
		currentValueTopN = 10;
	  }
    $.ajax({
		url: '/AJAXAlignment',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({
							'uniqID': uniqID,
							'kinase': kinase,
							'mutation': mutation,
							'results': results,
							'WS': currentValueWS,
							'topN': currentValueTopN,
						}),
		success: function (response){
            // document.getElementById('alignmentCard').innerHTML = response['filepath'];
			const image = document.getElementById('aliView');
			if (response['status'] != 'error') {
				image.setAttribute('src', response['filepath']);
				image.setAttribute('name', response['filepath'].split("/").pop());
				image.setAttribute('style', 'overflow: scroll;');
				image.setAttribute('alt', kinase+'/'+mutation);
				// alert (response['filepath']);
			}
			else {
				document.getElementById('alignmentButtonsCard').style.display = 'none';
				image.setAttribute('src', '');
				image.setAttribute('alt', 'Could not generate the alignment.');
				// image.setAttribute('style', 'height: 500px; width: 100%; object-fit: contain');
			}

		}
	});
}