/**
 * Menu
 */
 $("a.menu-icon").on("click", function(event) {
   var w = $(".menu");

   w.css({
     display: w.css("display") === "none"
      ? "block"
      : "none"
   });
 });

$("#homeimgcontainer")[0].innerHTML += "<div>서울, January 2019.</div>"

 $("script[type='math/tex']").replaceWith(
    function(){
      var tex = $(this).text();
      return katex.renderToString(tex, {displayMode: false});
  });
  
  $("script[type='math/tex; mode=display']").replaceWith(
    function(){
      var tex = $(this).text();
      return katex.renderToString(tex.replace(/%.*/g, ''), {displayMode: true});
  });
