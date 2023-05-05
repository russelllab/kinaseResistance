window.onload = function() {
    document.getElementById("formArea").style.display = "block";
    document.getElementById("slider-text").style.display = "none";
    document.getElementById("slider-box").style.display = "none";

    // Get the text area element
    var textarea = document.getElementById("inputMut");
    
    // Create a new AJAX request
    var xhr = new XMLHttpRequest();

    // Set the URL of the text file
    var url = "/static/examples/example1.txt";

    // Open the request
    xhr.open("GET", url);

    // Set the response type to "text"
    xhr.responseType = "text";

    // Handle the response
    xhr.onload = function() {
        if (xhr.status === 200) {
            // Set the text area's value to the response text
            textarea.value = xhr.responseText;
        } else {
            // Handle the error
            console.log("Error loading text file");
        }
    };

    // Send the request
    xhr.send();
    
    // Set the text area's value to the default text
    // textarea.value = defaultText;
};

function call_onclick() {
    document.getElementById("formArea").style.display = "none";
    document.getElementById("slider-text").style.display = "block";
    document.getElementById("slider-box").style.display = "block";

    var div = document.getElementById("headingElement");
    div.style.position = "absolute";
    div.style.top = "40%";
    div.style.left = "50%";
    div.style.transform = "translate(-40%, -50%) scale(1.5)";
}

function call_onclick_clear() {
    // alert(document.getElementById("inputMut").innerHTML)
    document.getElementById("inputMut").value = "";
}

function load_example(num) {
    var textarea = document.getElementById("inputMut");
    var xhr = new XMLHttpRequest();
    var url = "/static/examples/example" + num + '.txt';
    xhr.open("GET", url);
    xhr.responseType = "text";
    xhr.onload = function() {
        if (xhr.status === 200) {
            textarea.value = xhr.responseText;
        } else {
            console.log("Error loading text file");
        }
    };
    xhr.send();
}