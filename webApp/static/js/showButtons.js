function setActiveClass(dropdownId) {
    // Rename the button
    var button = document.getElementById("dropdownMenuButton");
    button.innerHTML = event.target.innerHTML;

    // Set active class to the clicked button
    var dropdown = document.getElementById(dropdownId);
    var links = dropdown.getElementsByTagName('a');

    for (var i = 0; i < links.length; i++) {
        links[i].addEventListener('click', function() {
        // Remove active class from all links
        for (var j = 0; j < links.length; j++) {
            links[j].classList.remove('active');
        }

        // Add active class to the clicked link
        this.classList.add('active');
        });
    }
}

function showButtons()
{
    $.ajax({
		url: '/AJAXButtons',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({'uniqID': 'ABCDEF',}),
		success: function (response){
            // document.getElementById('summaryCard').innerHTML = response['text'];
            var dropdown = document.getElementById("dropdownButton");
            var button = document.createElement("button");
            button.setAttribute("class", "btn btn-secondary btn-sm dropdown-toggle");
            button.setAttribute("type", "button");
            button.setAttribute("id", "dropdownMenuButton");
            button.setAttribute("data-toggle", "dropdown");
            button.setAttribute("aria-haspopup", "true");
            button.setAttribute("aria-expanded", "false");
            button.innerHTML = "Functional Information";
            dropdown.appendChild(button);

            var divElement = document.createElement("div");
            divElement.setAttribute("class", "dropdown-menu");
            divElement.setAttribute("aria-labelledby", "dropdownMenuButton");

            // buttonList = ['Sequence identity', 'Phosphorylation', 'Acetylation']
            buttonList = ['Functional Information','Sequence Identity','Activating',
                            'Deactivating','Activating & Deactivating','Resistance','Phosphorylation','Acetylation',
                            'Ubiquitination','Sumoylation','O-GlcNAc','Methylation']
            for (let i = 0; i < buttonList.length; i++){
                var buttonAtag = document.createElement("a");
                buttonAtag.setAttribute("class", "dropdown-item");
                // buttonAtag.setAttribute("href", "#");
                if (i == 0) {
                    buttonAtag.setAttribute("class", "dropdown-item active");
                }
                buttonAtag.onclick = function() {
                    setActiveClass("dropdownButton");
                };
                buttonAtag.innerHTML = buttonList[i];
                divElement.appendChild(buttonAtag);
            }
            dropdown.appendChild(divElement);

		}
	});
}