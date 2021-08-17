/*global google*/
/*global navigator*/

const GOOGLE_API_KEY = 'AIzaSyAdWz_SegPgJPe1nC6vgJ09ojSxB_KfaOQ';

var map = null;
var marker = null;

function get_location() {
  if (navigator.geolocation) {
    navigator.geolocation.getCurrentPosition(function(position) {
      var lat = position.coords.latitude;
      var lon = position.coords.longitude;
      var latLng = new google.maps.LatLng(lat, lon);
      
      if (marker)
        marker.setMap(null);
        
      marker = new google.maps.Marker({
        position: latLng,
        map: map,
      });
      
      //map.setCenter(latLng);
      display_sir()
    });
  }
}

function refit_sir() {
  if (!marker)
    return;

  var lat = marker.getPosition().lat();
  var lon = marker.getPosition().lng();
  var n_days = document.getElementById('n_days').value;
  var detection_rate = document.getElementById('detection_rate').value;
  var piq = document.getElementById('piq').value;
  var pai = document.getElementById('pai').value;
  var refit = "True"

  var sir_uri = `/api/sir/${lat}/${lon}/${n_days}/${detection_rate}/${piq}/${pai}/${refit}`;
  console.log(sir_uri)
  var plot_img = document.getElementById('plot_img');
  plot_img.setAttribute('src', sir_uri);
}

function display_sir() {
  if (!marker)
    return;

  var lat = marker.getPosition().lat();
  var lon = marker.getPosition().lng();
  var n_days = document.getElementById('n_days').value;
  var detection_rate = document.getElementById('detection_rate').value;
  var piq = document.getElementById('piq').value;
  var pai = document.getElementById('pai').value;
  var refit = "False"

  var sir_uri = `/api/sir/${lat}/${lon}/${n_days}/${detection_rate}/${piq}/${pai}/${refit}`;
  console.log(sir_uri)
  var plot_img = document.getElementById('plot_img');
  plot_img.setAttribute('src', sir_uri);
}

function search() {
    if(event.key == 'Enter') {
        display_sir();
    }
}    

function create_map() {
  map = new google.maps.Map(document.getElementById('map'), {
    center: new google.maps.LatLng(39.8283, -98.5795),
    zoom: 4,
    mapTypeId: 'roadmap',
    disableDoubleClickZoom: true,
    zoomControl: false,
    streetViewControl: false,
  });
  
  map.addListener('dblclick', function(e) {
    if (marker)
      marker.setMap(null);
    
    marker = new google.maps.Marker({
      position: e.latLng,
      map: map,
    });
    
    //map.setCenter(e.latLng)
    display_sir();
  });
  
  get_location();
}
