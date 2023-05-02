function makeModals(modal_title)
{
    // modal_title = 'results';
    // alert (modal_title);
    var modalID = document.getElementById('modalPosition');
    // alert(responseText);
    $.ajax({
		url: '/AJAXModals',
		type: 'post',
		dataType: 'json',
		data: JSON.stringify({'modalTitle': modal_title}),
		success: function (response){
            // alert (response['modalText']);
            // return 'hello'
            var div1 = document.createElement("div");
            div1.setAttribute("class", "modal fade");
            div1.setAttribute("id", modal_title+"Modal");
            div1.setAttribute("tabindex", "-1");
            var div2 = document.createElement("div");
            div2.setAttribute("class", "modal-dialog");
            var div3 = document.createElement("div");
            div3.setAttribute("class", "modal-content");
            var div4 = document.createElement("div");
            div4.setAttribute("class", "modal-header");
            var h5 = document.createElement("h5");
            h5.setAttribute("class", "modal-title");
            h5.setAttribute("id", modal_title+"ModalLabel");
            h5.innerHTML = modal_title;
            var button1 = document.createElement("button");
            button1.setAttribute("type", "button");
            button1.setAttribute("class", "close");
            button1.setAttribute("data-dismiss", "modal");
            button1.setAttribute("aria-label", "Close");
            var span1 = document.createElement("span");
            span1.setAttribute("aria-hidden", "true");
            span1.innerHTML = "&times;";
            button1.appendChild(span1);
            div4.appendChild(h5);
            div4.appendChild(button1);
            var div5 = document.createElement("div");
            div5.setAttribute("class", "modal-body");
            div5.setAttribute("id", modal_title+"ModalBody");
            // div5.innerHTML = "Modal body text goes here.";
            // div5.innerHTML = loadModalText(modal_title);
            // alert('hello'+loadModalText(modal_title));
            div5.innerHTML = response['modalText'];
            var div6 = document.createElement("div");
            div6.setAttribute("class", "modal-footer");
            var button2 = document.createElement("button");
            button2.setAttribute("type", "button");
            button2.setAttribute("class", "btn btn-secondary");
            button2.setAttribute("data-dismiss", "modal");
            button2.innerHTML = "Close";
            div6.appendChild(button2);
            div3.appendChild(div4);
            div3.appendChild(div5);
            div3.appendChild(div6);
            div2.appendChild(div3);
            div1.appendChild(div2);

            modalID.appendChild(div1);
        }
    });
}