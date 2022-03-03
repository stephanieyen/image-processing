var Student = {
  // please fill in your name and NetID
  //  your NetID is the part of your email before @princeton.edu
  name: "Stephanie Yen",
  netID: "sy7",
};

Student.updateHTML = function() {
  var studentInfo = this.name + " &lt;" + this.netID + "&gt;";
  document.getElementById("student").innerHTML = studentInfo;
};
