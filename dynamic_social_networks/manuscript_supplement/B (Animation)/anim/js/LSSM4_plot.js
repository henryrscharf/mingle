(function($) {
    $(document).ready(function() {
	
	$('#LSSM4_plot').scianimator({
	    'images': ['../anim/LSSM4_images/LSSM4_plot1.png', '../anim/LSSM4_images/LSSM4_plot2.png', '../anim/LSSM4_images/LSSM4_plot3.png', '../anim/LSSM4_images/LSSM4_plot4.png', '../anim/LSSM4_images/LSSM4_plot5.png', '../anim/LSSM4_images/LSSM4_plot6.png', '../anim/LSSM4_images/LSSM4_plot7.png', '../anim/LSSM4_images/LSSM4_plot8.png', '../anim/LSSM4_images/LSSM4_plot9.png', '../anim/LSSM4_images/LSSM4_plot10.png', '../anim/LSSM4_images/LSSM4_plot11.png', '../anim/LSSM4_images/LSSM4_plot12.png', '../anim/LSSM4_images/LSSM4_plot13.png', '../anim/LSSM4_images/LSSM4_plot14.png', '../anim/LSSM4_images/LSSM4_plot15.png', '../anim/LSSM4_images/LSSM4_plot16.png', '../anim/LSSM4_images/LSSM4_plot17.png', '../anim/LSSM4_images/LSSM4_plot18.png', '../anim/LSSM4_images/LSSM4_plot19.png', '../anim/LSSM4_images/LSSM4_plot20.png', '../anim/LSSM4_images/LSSM4_plot21.png', '../anim/LSSM4_images/LSSM4_plot22.png', '../anim/LSSM4_images/LSSM4_plot23.png', '../anim/LSSM4_images/LSSM4_plot24.png', '../anim/LSSM4_images/LSSM4_plot25.png', '../anim/LSSM4_images/LSSM4_plot26.png', '../anim/LSSM4_images/LSSM4_plot27.png', '../anim/LSSM4_images/LSSM4_plot28.png', '../anim/LSSM4_images/LSSM4_plot29.png', '../anim/LSSM4_images/LSSM4_plot30.png', '../anim/LSSM4_images/LSSM4_plot31.png', '../anim/LSSM4_images/LSSM4_plot32.png', '../anim/LSSM4_images/LSSM4_plot33.png', '../anim/LSSM4_images/LSSM4_plot34.png', '../anim/LSSM4_images/LSSM4_plot35.png', '../anim/LSSM4_images/LSSM4_plot36.png', '../anim/LSSM4_images/LSSM4_plot37.png', '../anim/LSSM4_images/LSSM4_plot38.png', '../anim/LSSM4_images/LSSM4_plot39.png', '../anim/LSSM4_images/LSSM4_plot40.png', '../anim/LSSM4_images/LSSM4_plot41.png', '../anim/LSSM4_images/LSSM4_plot42.png', '../anim/LSSM4_images/LSSM4_plot43.png', '../anim/LSSM4_images/LSSM4_plot44.png', '../anim/LSSM4_images/LSSM4_plot45.png', '../anim/LSSM4_images/LSSM4_plot46.png', '../anim/LSSM4_images/LSSM4_plot47.png', '../anim/LSSM4_images/LSSM4_plot48.png', '../anim/LSSM4_images/LSSM4_plot49.png', '../anim/LSSM4_images/LSSM4_plot50.png', '../anim/LSSM4_images/LSSM4_plot51.png', '../anim/LSSM4_images/LSSM4_plot52.png', '../anim/LSSM4_images/LSSM4_plot53.png', '../anim/LSSM4_images/LSSM4_plot54.png', '../anim/LSSM4_images/LSSM4_plot55.png', '../anim/LSSM4_images/LSSM4_plot56.png', '../anim/LSSM4_images/LSSM4_plot57.png', '../anim/LSSM4_images/LSSM4_plot58.png', '../anim/LSSM4_images/LSSM4_plot59.png', '../anim/LSSM4_images/LSSM4_plot60.png', '../anim/LSSM4_images/LSSM4_plot61.png', '../anim/LSSM4_images/LSSM4_plot62.png', '../anim/LSSM4_images/LSSM4_plot63.png', '../anim/LSSM4_images/LSSM4_plot64.png', '../anim/LSSM4_images/LSSM4_plot65.png', '../anim/LSSM4_images/LSSM4_plot66.png', '../anim/LSSM4_images/LSSM4_plot67.png', '../anim/LSSM4_images/LSSM4_plot68.png', '../anim/LSSM4_images/LSSM4_plot69.png', '../anim/LSSM4_images/LSSM4_plot70.png', '../anim/LSSM4_images/LSSM4_plot71.png', '../anim/LSSM4_images/LSSM4_plot72.png', '../anim/LSSM4_images/LSSM4_plot73.png', '../anim/LSSM4_images/LSSM4_plot74.png', '../anim/LSSM4_images/LSSM4_plot75.png', '../anim/LSSM4_images/LSSM4_plot76.png', '../anim/LSSM4_images/LSSM4_plot77.png', '../anim/LSSM4_images/LSSM4_plot78.png', '../anim/LSSM4_images/LSSM4_plot79.png', '../anim/LSSM4_images/LSSM4_plot80.png', '../anim/LSSM4_images/LSSM4_plot81.png', '../anim/LSSM4_images/LSSM4_plot82.png', '../anim/LSSM4_images/LSSM4_plot83.png', '../anim/LSSM4_images/LSSM4_plot84.png', '../anim/LSSM4_images/LSSM4_plot85.png', '../anim/LSSM4_images/LSSM4_plot86.png', '../anim/LSSM4_images/LSSM4_plot87.png', '../anim/LSSM4_images/LSSM4_plot88.png', '../anim/LSSM4_images/LSSM4_plot89.png', '../anim/LSSM4_images/LSSM4_plot90.png', '../anim/LSSM4_images/LSSM4_plot91.png', '../anim/LSSM4_images/LSSM4_plot92.png', '../anim/LSSM4_images/LSSM4_plot93.png', '../anim/LSSM4_images/LSSM4_plot94.png', '../anim/LSSM4_images/LSSM4_plot95.png', '../anim/LSSM4_images/LSSM4_plot96.png', '../anim/LSSM4_images/LSSM4_plot97.png'],
	    'width': 650,
	    'delay': 160,
	    'loopMode': 'loop',
 'controls': ['first', 'previous', 'play', 'next', 'last', 'loop', 'speed']
	});
	$('#LSSM4_plot').scianimator('play');
    });
})(jQuery);