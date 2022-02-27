package fractureW20;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Color3f;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import com.jogamp.opengl.util.gl2.GLUT;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.FileSelect;
import mintools.swing.HorizontalFlowPanel;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.EasyViewer;
import mintools.viewer.FlatMatrix4d;
import mintools.viewer.Interactor;
import mintools.viewer.SceneGraphNode;

/**
 * Fracture and continuum elastic simulation application
 * @author kry
 */
public class FractureApp implements SceneGraphNode, Interactor {

    private EasyViewer ev;
    
    public FEMSystem system;

    /**
     * Entry point for application
     * @param args
     */
    public static void main(String[] args) {
        new FractureApp();        
    }
        
    /**
     * Creates the application / scene instance
     */
    public FractureApp() {
    	T.getBackingMatrix().setIdentity();
        system = new FEMSystem();
        ev = new EasyViewer( "FEM + Fracture", this, new Dimension(640,480), new Dimension(640,480) );
        // we add ourselves as an interactor to set up mouse and keyboard controls
        ev.addInteractor(this);
    }

    /** Transform from world to screen */
    private FlatMatrix4d T = new FlatMatrix4d();

    /**
     * Converts from screen coordinates to world coordinates
     * @param x screen x
     * @param y screen y
     * @param p image coordinates
     */
    private void setPoint( int x, int y, Point2d p ) {
        Point3d tmp = new Point3d(x,y,0);
        Matrix4d M = new Matrix4d();
        M.invert(T.getBackingMatrix());
        M.transform( tmp );
        p.set( tmp.x, tmp.y );   
    }
    
    @Override
    public void init(GLAutoDrawable drawable) {
        GL gl = drawable.getGL();
        gl.glEnable( GL.GL_BLEND );
        gl.glBlendFunc( GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA );
        gl.glEnable( GL.GL_LINE_SMOOTH );
        gl.glEnable( GL2.GL_POINT_SMOOTH );
        system.init(drawable);
    }
        
    @Override
    public void display(GLAutoDrawable drawable) {

        GL2 gl = drawable.getGL().getGL2();        
        gl.glClearColor( 1, 1, 1, 1);
        gl.glClear( GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT );
        
        // First advance the system (if it is running or wants to be stepped)
        if ( run.getValue() || stepRequest ) {
            for ( int i = 0; i < substeps.getValue(); i++ ) {
                if ( ! system.updateParticles( stepsize.getValue() / substeps.getValue() )) {
                    run.setValue( false );
                    break;
                }
            }
        }
                
        gl.glDisable( GL2.GL_LIGHTING );

        // We're doing 2D drawing only, so we'll 
        // call this helper function to set up a 2D projection 
        // where the origin is at the top left corner and the 
        // units are pixels
        EasyViewer.beginOverlay(drawable);
        
        gl.glPushMatrix();
                
	        windowWidth = drawable.getSurfaceWidth();
	        windowHeight = drawable.getSurfaceHeight();
	        double sw = scale.getValue() * windowWidth / worldWidth;
	        double sh = scale.getValue() * windowHeight / worldHeight;
	        double s = Math.min(sw,sh);        
	        double x = posx.getValue(); // extra translation (world center)
	        double y = posy.getValue();                
	        // create this as a transform to use for mouse picking and interaction        
	        Matrix4d M = T.getBackingMatrix();
	        M.setIdentity();
	        M.m00 = s;
	        M.m11 = s;        
	        M.m03 = x + windowWidth /2 - s * worldCenterx; 
	        M.m13 = y + windowHeight/2 - s * worldCentery; 
	                
	        gl.glMultMatrixd( T.asArray(),0 );
	           
	        
	        system.display( drawable );
        
        
        gl.glPopMatrix();
        
        // Finally we'll display a string with useful information
        // about the system and the current stepping      
        if ( system.showCommentsAndParameters.getValue() ) {
	        String text = system.toString();
	        if ( displayStepInformation.getValue() ) {
	        	text += "\n" + "h = " + stepsize.getValue() + "\n" +
	        			"substeps = " + (int)(double) substeps.getValue();
	        }
	        gl.glColor3f(0.5f,0.5f,0.5f);
	        EasyViewer.printTextLines( drawable, text, 10, 15, 18, GLUT.BITMAP_9_BY_15 );
        }
        EasyViewer.endOverlay(drawable);    

        // If we're recording, we'll save the step to an image file.
        // we'll also clear the step request here.
        if ( run.getValue() || stepRequest ) {
            stepRequest = false;        
            if ( record.getValue() ) {
                // write the frame
                File file = new File( "stills/" + dumpName + format.format(nextFrameNum) + ".png" );                                             
                nextFrameNum++;
                file = new File(file.getAbsolutePath().trim());
                ev.snapshot(drawable, file);
            }
        }
    }
    
    /**
     * Computes a mass based on the intensity, darker is heavier.  Uses standard colour to intensity conversion weights.
     * @return a value between 0 and 1
     */
    public double getColourMass( Color3f c ) {
        return 1 - (0.3 * c.x + 0.59 * c.y + 0.11 * c.z);
    }
    
    private void loadSystemFromImage( ) {
    	File f = FileSelect.select("png", "image", "load", "./dataFracture", true);
    	if ( f != null ) {
    		loadSystemFromImage( f.getAbsolutePath() );
    	}
    }
    
    /**
     * Loads the specified image, clearing the old system, and resets viewing parameters.
     * @param filename
     */
    private void loadSystemFromImage( String filename ) {
        systemClear();
        system.name = filename;
        // scale up each pixel by factor of 10 to try to get these systems to have similar 
        // behaviour to the systems generated with Triangle
        ImageBlocker blocker = new ImageBlocker( new File(filename), 0.05f, 10, randomPointDisplacement.getValue(), fourTrianglesPerPixel.getValue() );
		system.femSprings.addAll( blocker.getTriangles() );
		system.particles.addAll( blocker.getParticles() );
		setViewingPositionAndScale();
		for ( Particle p : system.particles ) {
			p.mass = 0;
		}
		for ( FEMTriangle t : system.femSprings ) {
			// let the darkness of the colour be the density
			// TODO: note the arbitrary mass multiplier here
			double m = 0.001 * getColourMass( t.colour ) * t.area / 3;
			t.A.mass += m;
			t.B.mass += m;
			t.C.mass += m;
		}
		system.identifyBoundaries();
    }

    private void loadSystemFromTriangles() {
	    File f = FileSelect.select("ele", "triangles", "load", "./dataFracture", true);                	
		if ( f != null ) {
			String name = f.getAbsolutePath();
			loadSystemFromTriangles( name );
		}
    }
    
    /** 
     * Loads from triangle file 
     */
    private void loadSystemFromTriangles( String filename ) {
        systemClear();
        system.name = filename;
		// Can use M to scale and translate wahtever is being loaded as a triangle mesh
		Matrix3d M = new Matrix3d();
		M.setIdentity();
		String filenameroot = filename.substring(0, filename.length()-4 );
		Triangle.loadTriangles( system, filenameroot, M );
		setViewingPositionAndScale();
		system.identifyBoundaries();
    }
    
    /**
     * Sets the viewing position and scale based on positions of particles in the system.
     */
    public void setViewingPositionAndScale() {
    	Point2d min = new Point2d( system.particles.get(0).p );
		Point2d max = new Point2d( system.particles.get(0).p );
		for ( Particle p : system.particles ) {
			if ( p.p.x < min.x ) min.x = p.p.x;
			if ( p.p.y < min.y ) min.y = p.p.y;
			if ( p.p.x > max.x ) max.x = p.p.x;
			if ( p.p.y > max.y ) max.y = p.p.y;
		}
		worldWidth = max.x - min.x;
		worldHeight = max.y - min.y;
		worldCenterx = 0.5*(min.x + max.x);
		worldCentery = 0.5*(min.y + max.y);
		posx.setValue(0);
		posy.setValue(0);
    }
    
    
    /** 
     * boolean to signal that the system was stepped and that a 
     * frame should be recorded if recording is enabled
     */
    private boolean stepRequest = false;
        
    /**
     * Base name of images to save
     */
    private String dumpName = "img";
    
    /**
     * Index for the frame number we're saving.
     */
    private int nextFrameNum = 0;
    
    /**
     * For formating the image file name when recording frames
     */
    private NumberFormat format = new DecimalFormat("00000");
        
    private BooleanParameter record = new BooleanParameter( "record each step to image file (press ENTER in canvas to toggle)", false );

    private BooleanParameter run = new BooleanParameter( "simulate (press SPACE in canvas to toggle)", false );
    
    private DoubleParameter stepsize = new DoubleParameter( "step size", 0.01, 1e-5, 1 );
    
    private DoubleParameter substeps = new DoubleParameter( "sub steps (integer)", 1, 1, 100);
        
    private BooleanParameter displayStepInformation = new BooleanParameter( "display stepping information overlay", true );
    
    private DoubleParameter randomPointDisplacement = new DoubleParameter( "random point displacement for image based load", 0.1, 0, 0.2 );
    private BooleanParameter fourTrianglesPerPixel = new BooleanParameter( "4 (otherwise 2) tri per pixel", true );
    
    // parameters and variables for for scaling and translating the window
    private DoubleParameter scale = new DoubleParameter("scale scene",.9, 0.1, 10);
    private DoubleParameter posx = new DoubleParameter("x translation", 0, -10000, 10000 );
    private DoubleParameter posy = new DoubleParameter("y translation", 0, -10000, 10000 );
    private double worldWidth = 1.0;
    private double worldHeight = 1.0;
    private double worldCenterx = 0.5;
    private double worldCentery = 0.5;
    private int windowWidth;
    private int windowHeight;
    
    private Point prevMousePoint = new Point();    
    private Point2d mousePoint = new Point2d();   
    
    @Override
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        
        vfp.add( record.getControls() );  
        
        JButton loadTrianglesELE = new JButton("Load from .ele file");
        loadTrianglesELE.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				loadSystemFromTriangles();
			}
		});
        vfp.add( loadTrianglesELE );
        JButton loadTrianglesPNG = new JButton("Load from .png file");
        loadTrianglesPNG.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				loadSystemFromImage();
			}
		});
        vfp.add( loadTrianglesPNG );
        vfp.add( randomPointDisplacement.getSliderControls(false) );
        vfp.add( fourTrianglesPerPixel.getControls() );
        
        HorizontalFlowPanel hfp2 = new HorizontalFlowPanel();
        hfp2.add( record.getControls() );
        JButton res1 = new JButton("640x360");
        hfp2.add( res1);
        res1.addActionListener( new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ev.glCanvas.setSize( 640, 360 );
                ev.frame.setSize( ev.frame.getPreferredSize() );
            }
        });        
        JButton res2 = new JButton("1280x720");
        hfp2.add( res2);
        res2.addActionListener( new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {                
                ev.glCanvas.setSize( 1280, 720 );
                ev.frame.setSize( ev.frame.getPreferredSize() );

            }
        });                
        vfp.add( hfp2.getPanel() );
        
        VerticalFlowPanel vfp0 = new VerticalFlowPanel();
        vfp0.setBorder( new TitledBorder("Numerical Integration Controls"));
        vfp0.add( run.getControls() );        
        vfp0.add( stepsize.getSliderControls(true) );
        vfp0.add( substeps.getControls() );
        vfp0.add( displayStepInformation.getControls() );
        CollapsiblePanel cp0 = new CollapsiblePanel( vfp0.getPanel() );
        //cp0.collapse();
        vfp.add( cp0 );
        
        vfp.add( system.getControls() );
                
        VerticalFlowPanel vfpv = new VerticalFlowPanel();
        vfpv.setBorder( new TitledBorder("window content scaling controls") );
        vfpv.add( scale.getSliderControls(true) );
        vfpv.add( posx.getSliderControls(false) );
        vfpv.add( posy.getSliderControls(false) );
        CollapsiblePanel vcp = new CollapsiblePanel(vfpv.getPanel());
        vcp.collapse();
        vfp.add( vcp );
        
        return vfp.getPanel();
    }
        
    private void systemReset() {
    	if ( system.name.endsWith("png") ) {
    		loadSystemFromImage( system.name );
    	} else if ( system.name.endsWith("ele") ){
    		loadSystemFromTriangles( system.name );    		
    	}
    }
    
    /**
     * Clears the system and resets viewing
     */
    private void systemClear() {
        posx.setValue(0.0);
        posy.setValue(0.0);
        scale.setValue(0.9);
        system.clear();
    }
    
    @Override
    public void attach(Component component) {
        component.addMouseMotionListener( new MouseMotionListener() {
            @Override
            public void mouseDragged(MouseEvent e) {            	
            	setPoint( e.getPoint().x, e.getPoint().y, mousePoint );
            	// update mouse spring
            	system.mouseSpring.mousePos.set( mousePoint );
                // handle translation in interface
                if ( (e.getModifiers() & InputEvent.BUTTON2_MASK) != 0) { // button3 ) {
                    posx.setValue( posx.getValue() + e.getPoint().x - prevMousePoint.x );
                    posy.setValue( posy.getValue() + e.getPoint().y - prevMousePoint.y );
                }
                // handle scaling in interface
                if ( (e.getModifiers() & InputEvent.BUTTON3_MASK) != 0 ) {
                    scale.setValue( scale.getValue() * Math.pow(1.002, e.getPoint().y - prevMousePoint.y ));
                }
                prevMousePoint.setLocation( e.getPoint() );
            }
            @Override
            public void mouseMoved(MouseEvent e) {
                // do nothing
            }
        } );
        component.addMouseListener( new MouseListener() {
            @Override
            public void mouseClicked(MouseEvent e) {
            	double s = scale.getValue() * windowWidth / worldWidth;
            	// Close should be determined by the current screen scaling
            	// We'll toggle anything with 10 pixels.
            	setPoint( e.getPoint().x, e.getPoint().y, mousePoint );
            	for ( Particle p : system.particles ) {
            		if ( p.p.distance( mousePoint ) < 10/s ) {
            			p.pinned = ! p.pinned;
            		}
            	}
            }
            @Override
            public void mouseEntered(MouseEvent e) {
            	// do nothing
            }
            @Override
            public void mouseExited(MouseEvent e) {
            	// do nothing
            }
            @Override
            public void mousePressed(MouseEvent e) {
            	prevMousePoint.setLocation( e.getPoint() );
                setPoint( e.getPoint().x, e.getPoint().y, mousePoint ); 
            	system.mouseSpring.triangle = null;
                if ( (e.getModifiers() & InputEvent.BUTTON1_MASK) != 0 ) {
                	for ( FEMTriangle triangle : system.femSprings ) {
                		if ( triangle.isInside( mousePoint, 0, system.mouseSpring.barycentricCoords ) ) {
                			system.mouseSpring.triangle = triangle;
                			system.mouseSpring.mousePos.set( mousePoint );
                			system.useMouseSpring = true;
                			break;
                		}
                	}
                }
            }
            @Override
            public void mouseReleased(MouseEvent e) {
            	system.useMouseSpring = false;
            }
        } );
        component.addKeyListener( new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                if ( e.getKeyCode() == KeyEvent.VK_SPACE ) {
                    run.setValue( ! run.getValue() ); 
                } else if ( e.getKeyCode() == KeyEvent.VK_S ) {                    
                    stepRequest = true;
                } else if ( e.getKeyCode() == KeyEvent.VK_R ) {
                    systemReset();
                    //system.resetParticles();
                } else if ( e.getKeyCode() == KeyEvent.VK_C ) {                   
                    systemClear();
                } else if ( e.getKeyCode() == KeyEvent.VK_1 ) {
                	loadSystemFromTriangles();                	
                } else if ( e.getKeyCode() == KeyEvent.VK_2 ) {
                	loadSystemFromImage();
                } else if ( e.getKeyCode() == KeyEvent.VK_ESCAPE ) {
                    // quit the program
                    ev.stop();
                } else if ( e.getKeyCode() == KeyEvent.VK_ENTER ) {
                    // toggle recording of steps to png files
                    record.setValue( ! record.getValue() );
                }
                if ( e.getKeyCode() != KeyEvent.VK_ESCAPE ) ev.redisplay();
            }
        } );
    }
    
}
