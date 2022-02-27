package fractureW20;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.vecmath.Color3f;
import javax.vecmath.Point2d;
import javax.vecmath.Tuple2d;
import javax.vecmath.Vector2d;

/**
 * Creates rigid bodies from an image
 * @author kry
 */
public class ImageBlocker {

    /** image data in format ARGB */
    private int[] imageData;
        
    /** image width */
    int width;
    
    /** image height */
    int height;
    
    private ArrayList<FEMTriangle> triangles = new ArrayList<FEMTriangle>();

    /**
     * Random numbers for perturbing point positions slightly
     */
    private Random r = new Random();

    private void randomDisplacement( double size, Tuple2d u ) {
    	u.x = r.nextDouble() * size - size/2;
    	u.y = r.nextDouble() * size - size/2;
    }
    
    /**
     * Creates a set of rigid bodies form the given image
     * @param filename
     * @param epsilon  threshold for detecting white, zero means white must be white
	 * @parame scale   for scaling up the dimensions to help normalize the model for default simulation parameters
	 * @parame randomDisplacementScale  can help break symmetry problems in CD with a bit of random displacement of DOFs
	 * @params fourTrianglesPerPixel  4 if true, 2 triangles if false with random selection of which
     */
    public ImageBlocker( File file, float epsilon, double scale, double randomDisplacementScale, boolean fourTrianglesPerPixel ) {
        this.epsilon = epsilon;
        // put a bit of noise on the points to try to help slightly buggy collision detection
        // avoid so many perfectly aligned cases.  But let's make it deterministic.
        try {
            BufferedImage img = ImageIO.read( file );
            width = img.getWidth();
            height = img.getHeight();
            imageData = new int[width*height];
            img.getRGB( 0, 0, width, height, imageData, 0, width );
            Color3f colour = new Color3f();
            for ( int x = 0; x < width; x++ ) {
                for ( int y = 0; y < height; y++ ) {
                    getColour( colour, x, y );
                    if ( isWhite( colour ) ) continue;
                    boolean isBlue = (colour.x == colour.y && colour.x < colour.z); 
                    if ( fourTrianglesPerPixel ) {
	                    // make 4 triangles with node at center.
	                    // find boundary edges later!
	
	                    // A D
	                    //  E
	                    // B C
	
	                    Particle A = getParticle( new Point2d( x, y ) );
	                    Particle B = getParticle( new Point2d( x, y+1 ) );
	                    Particle C = getParticle( new Point2d( x+1, y+1 ) );
	                    Particle D = getParticle( new Point2d( x+1, y ) );
	                    Particle E = getParticle( new Point2d( x+0.5, y+0.5 ) );
	                    if ( isBlue ) {
	                    	A.pinned = true;
	                    	B.pinned = true;
	                    	C.pinned = true;
	                    	D.pinned = true;
	                    	E.pinned = true;
	                    }
	                    // might seem CW instead of CCW, but image y is different than scene y
	                    triangles.add( new FEMTriangle( A, D, E, colour ) );
	                    triangles.add( new FEMTriangle( A, E, B, colour ) );
	                    triangles.add( new FEMTriangle( B, E, C, colour ) );
	                    triangles.add( new FEMTriangle( D, C, E, colour ) );     
                    } else {
	                    Particle A = getParticle( new Point2d( x, y ) );
	                    Particle B = getParticle( new Point2d( x, y+1 ) );
	                    Particle C = getParticle( new Point2d( x+1, y+1 ) );
	                    Particle D = getParticle( new Point2d( x+1, y ) );
	                    if ( isBlue ) {
	                    	A.pinned = true;
	                    	B.pinned = true;
	                    	C.pinned = true;
	                    	D.pinned = true;
	                    }
	                    // might seem CW instead of CCW, but image y is different than scene y
	                    if ( r.nextBoolean() ) {
		                    triangles.add( new FEMTriangle( A, D, B, colour ) );
		                    triangles.add( new FEMTriangle( D, C, B, colour ) );
	                    } else {
		                    triangles.add( new FEMTriangle( A, D, C, colour ) );
		                    triangles.add( new FEMTriangle( B, A, C, colour ) );   
	                    }
                    }
                }
            }
            
            r.setSeed(0); 
            Vector2d u = new Vector2d();
            for ( Particle p : point2ParticleMap.values() ) {
            	randomDisplacement( randomDisplacementScale, u );
            	p.p0.add( u );
            	p.p0.scale( scale );
            	p.p.set(p.p0);
            }
            for ( FEMTriangle t : triangles ) {
            	t.computeAreaAndDm();
            }
            
        } catch ( Exception e ) {
            System.err.println("Problems processing image " + file );
            e.printStackTrace();
        }
    }
    
    public Collection<Particle> getParticles() {
    	return point2ParticleMap.values();
    }
    
    public Collection<FEMTriangle> getTriangles() {
    	return triangles;
    }
    
    private HashMap<Point2d,Particle> point2ParticleMap = new HashMap<Point2d,Particle>();
    
    private Particle getParticle( Point2d p ) {
    	if ( point2ParticleMap.containsKey(p) ) {
    		return point2ParticleMap.get(p);
    	}
    	Particle A = new Particle( p.x, p.y, 0, 0 );
    	point2ParticleMap.put( p, A );
    	return A;
    }
    
    /** 
     * Gets the colour of the specified pixel
     * @param x
     * @param y
     */
    private void getColour( Color3f colour, int x, int y ) {
        int data = imageData[y*width+x];
        colour.x = ((data >> 16) & 0x0ff) / 255.0f;
        colour.y = ((data >> 8) & 0x0ff) / 255.0f;
        colour.z = ((data >> 0) & 0x0ff) / 255.0f;        
    }

    /** epsilon for checking for white pixels */
    private float epsilon = 0;
    
    /**
     * @param colour
     * @return true if the colour provided is white
     */
    private boolean isWhite( Color3f colour ) {
        final Color3f white = new Color3f(1,1,1);
        return colour.epsilonEquals(white, epsilon );
    }
    
}
