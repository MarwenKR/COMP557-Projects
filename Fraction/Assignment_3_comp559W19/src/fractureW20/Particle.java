package fractureW20;

import java.util.ArrayList;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

/**
 * Particle, i.e., Degree of Freedom, or DOF for the finite element model
 * @author kry
 */
public class Particle {
    
    /** true means that the particle can not move */
	public boolean pinned = false;
        
    /** The mass of the particle */
	public double mass = 1;
    
    /** current position of the particle */
    public Point2d p = new Point2d();
    
    /** current velocity of the particle */
    public Vector2d v = new Vector2d();
    
    /** initial position of the particle */
    public Point2d p0 = new Point2d();
    
    /** initial velocity of the particle */
    public Vector2d v0 = new Vector2d();
    
    /** force acting on this particle */
    public Vector2d f = new Vector2d();
    
    /** position differential for implicit integration*/
    public Vector2d dx = new Vector2d();
    
    /** force differential for implicit integration */
    public Vector2d df = new Vector2d();
    
    /**df/dv for implicit collision */
    public Vector2d dfdv_col = new Vector2d();

    /** temporary velocity of the particle */
    Vector2d vTmp = new Vector2d();
    
    /**
     * A list of adjacent triangles, necessary for processing fracture
     */
    public ArrayList<FEMTriangle> tris = new ArrayList<FEMTriangle>();
    
    public ArrayList<Vector2d> fplus = new ArrayList<Vector2d>();
    
    public ArrayList<Vector2d> fminus = new ArrayList<Vector2d>();
    
    /** The separtaion tensor for deciding fracture points and directions */
    Matrix2d separationTensor = new Matrix2d();
    
    public int addTriangle( FEMTriangle tri ) {
    	int i = tris.size();
    	tris.add( tri );    	
    	fplus.add( new Vector2d() );
    	fminus.add( new Vector2d() );
    	return i;
    }
    
    /** 
     * Compute the separation tensor.
     */
    public void computeSeparationTensor() {
    	separationTensor.zero();
    	
    	// TODO: Objective 3: Compute separation tensor using fPlus and fMinus lists
    
    	for (Vector2d f_plus: fplus) 
    		separationTensor.rank1(f_plus.length()!=0.0 ? 1.0/f_plus.length():0.0, f_plus);
    	    	
    	for (Vector2d f_minus: fminus) 	
    		separationTensor.rank1(f_minus.length() !=0.0 ? -1.0/f_minus.length():0.0, f_minus);
    	
    	// remove unbalanced tensile load
    	Vector2d F_plus = new Vector2d();
		F_plus.set(0.0,0.0);
    	for (Vector2d f_plus: fplus)
    		F_plus.add(f_plus);
    	
		separationTensor.rank1(F_plus.length()!=0.0 ? -1.0/F_plus.length():0.0, F_plus);
    	
		// remove unbalanced compressive load
		Vector2d F_minus = new Vector2d();
		F_minus.set(0.0,0.0);
    	for (Vector2d f_minus: fminus)
    		F_minus.add(f_minus);	
		separationTensor.rank1(F_minus.length()!=0.0 ? 1.0/F_minus.length():0.0, F_minus);
		
		//multiply by 1/2
		separationTensor.scale(0.5);
		
		separationTensor.evd();

    }
    
    public void drawSeparationTensor( GLAutoDrawable drawable, double s ) {
    	GL2 gl = drawable.getGL().getGL2();
    	gl.glBegin( GL.GL_LINES );
		double s1 = s * separationTensor.ev1;
		double s2 = s * separationTensor.ev2;
		if ( separationTensor.ev1 < 0 ) {
			gl.glColor3f(0.75f,0,0);
		} else {
			gl.glColor3f(0,0.75f,0);
		}
		gl.glVertex2d( p.x - separationTensor.v1.x * s1, p.y - separationTensor.v1.y * s1 );
		gl.glVertex2d( p.x + separationTensor.v1.x * s1, p.y + separationTensor.v1.y * s1 );
		if ( separationTensor.ev2 < 0 ) {
			gl.glColor3f(0.75f,0,0);
		} else {
			gl.glColor3f(0,0.75f,0);
		}
		gl.glVertex2d( p.x - separationTensor.v2.x * s2, p.y - separationTensor.v2.y * s2 );
		gl.glVertex2d( p.x + separationTensor.v2.x * s2, p.y + separationTensor.v2.y * s2 );

    	gl.glEnd();
    }
    
    /**
     * Creates a partilce copy, useful for fracture events
     * @param p
     */
    public Particle( Particle A ) {
    	p.set( A.p );
    	p0.set( A.p0 );
    	v.set( A.v );
    	v0.set( A.v0 );
    	vTmp.set( v0 ); // used in repulsion? not exactly sure what this is for
    	// Note that mass is computed elsewhere
    }
    
    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    public Particle( double x, double y, double vx, double vy ) {
        p0.set(x,y);
        v0.set(vx,vy);
        reset();
    }
    
    /**
     * Resets the position of this particle
     */
    public void reset() {
        p.set(p0);
        v.set(v0);
        f.set(0,0);
        vTmp.set(v0);
    }
    
    /**
     * Adds the given force to this particle.
     * Note that you probably want to set the force to zero 
     * before accumulating forces. 
     * @param force
     */
    public void addForce( Vector2d force ) {
        f.add(force);
    }
    
    /**
     * Computes the distance of a point to this particle
     * @param x
     * @param y
     * @return the distance
     */
    public double distance( double x, double y ) {
        Point2d tmp = new Point2d( x, y );
        return tmp.distance(p);
    }
    
    public int index;
   
}
