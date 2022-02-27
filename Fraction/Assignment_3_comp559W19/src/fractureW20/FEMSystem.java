package fractureW20;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.HorizontalFlowPanel;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.SceneGraphNode;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

/**
 * Implementation of a simple particle system.
 * 
 * Note that the particle system implements Function, that is, it evaluates
 * its derivatives to return to the step method which is called by implementations
 * of the Integrator interface.
 * 
 * Note also that it is actually the updateParticles method in this class which 
 * should be calling the Integrator step method! 
 * 
 * @author kry
 */
public class FEMSystem implements SceneGraphNode, Filter, MatrixMult {
    
    /** the particle list */
    public List<Particle> particles = new LinkedList<Particle>();
    
    /** particles for applying traction (i.e., bottom row of block tests) */
    public List<Particle> tractionParticles = new LinkedList<Particle>();

    /** the spring list for collisions */
    public List<Edge> collidableEdges = new LinkedList<Edge>();
    
    /** leaf springs connect 3 particles with a hinge like spring */
    public List<FEMTriangle> femSprings = new LinkedList<FEMTriangle>();
    
    /** List of collision objects **/
    public List<Collision> collisionList = new LinkedList<Collision>();
      
    /** The name of the particle system. */
    public String name = "";
    

    MouseSpring mouseSpring;
    boolean useMouseSpring = false;
    
    /**
     * Creates an empty particle system
     */
    public FEMSystem() {
    	mouseSpring = new MouseSpring();
    }
    
    /**
     * Resets the positions of all particles to their initial states.
     * Probably don't want to use this after fracture events as the
     * original fracture state will not be restored
     */
    public void resetParticles() {
    	collisionList.clear();
        for ( Particle p : particles ) {
            p.reset();
        }
        for ( FEMTriangle tri : femSprings ) {
        	tri.reset(); // try to get rid of NaNs
        }
        time = 0;
    }
    
    /**
     * Deletes all particles, and as such removes all springs too.
     */
    public void clear() {      
    	collisionList.clear();
        particles.clear();
        tractionParticles.clear();
        collidableEdges.clear();
        femSprings.clear();
        name = "";
    }    
    
    public double time = 0;
    
    private double stepSize;
    

    /**
     * Identifies the boundary edges of the current set of FEM triangles.
     * The method loops over all triangles adding edges to a set, and then 
     * removing those edges if the edge is found again in an adjacent triangle. 
     * This is needed on loading a new system.  You'll want to do something 
     * more efficient on fixing up boundaries after fracture events.
     */
    void identifyBoundaries() {
    	collidableEdges.clear();
		final HashSet<Edge> boundaries = new HashSet<Edge>();
		boundaries.clear();				
		for ( FEMTriangle t : femSprings ) {
			Edge e1 = new Edge( t.A, t.B, t );
			if ( ! boundaries.remove( e1 ) ) boundaries.add( e1 );
			Edge e2 = new Edge( t.B, t.C, t );
			if ( ! boundaries.remove( e2 ) ) boundaries.add( e2 );
			Edge e3 = new Edge( t.C, t.A, t );
			if ( ! boundaries.remove( e3 ) ) boundaries.add( e3 );
		}
		collidableEdges.addAll( boundaries ); 
    }
    
    
    
    
    /**
     * Advances time and updates the position of all particles
     * @param h 
     * @return true if update was successful
     */
    public boolean updateParticles( double h ) {
        boolean resultOK = true;
        Vector2d tmp = new Vector2d();
        stepSize = h;

        // set the global element properties
        double E = YoungModulus.getValue();
        double nu = PoissonRatio.getValue();
        FEMTriangle.mu = E/2/(1+nu);
        FEMTriangle.lambda = E*nu/(1+nu)/(1-2*nu);
        mouseSpring.k = mouseSpringk.getValue();
        
        collisionList.clear();
        
        // ensure proper indexing for global deltaV vector
        for (int index = 0; index < particles.size(); index++) {
        	particles.get(index).index = index;
        }
          
        if ( doCollisionDetection.getValue() ) {
	        // process collision
        	// Note that this would be better as a loop over all boundary vertices
        	for (int i = 0; i < collidableEdges.size(); i++) {
        		FEMTriangle t1 = collidableEdges.get(i).triangle;
	        	for (int j = i + 1; j < collidableEdges.size(); j++) {
	        		FEMTriangle t2 = collidableEdges.get(j).triangle;
	        		double eps = 1e-10;
	        		Collision c = Collision.collisionDetect( t1, t2, eps );
	        		if ( c != null ) {
	        			c.coefficent_viscous = collisionViscousCoefficient.getValue();
	        			c.relNormalVelThres = relNormalVelThresh.getValue();
	        			c.coefficent_repulsive = collisionRepulsion.getValue(); 
	        			c.updateCollisionResponse();
	        			collisionList.add( c );
	        		}	        		
	        	}
	        }        	
        }
               		
        // Update the velocity of the particles as per symplectic Euler
        if ( ! implicit.getValue() ) {
        	
        	computeForces();
	        for ( Particle p : particles ) {
	            if ( p.pinned ) {            
	                p.f.set(0,0); // just to make sure!
	                p.v.set(0,0);
	            } else {
	                tmp.scale( h / p.mass, p.f );
	                p.v.add( tmp );            
	            }
	        }
        } else {
        	
        	// TODO: Objective 2: advance the state of the system with implicit integration
        	// it is fine to do only 1 linear solve rather than Newton iterations (i.e.,
        	// you can ignore the newtonIterations parameter)
        	
        	init(); // initialize useful members you want to use here based on the number of DOFs
        	
        	// get velocity
        	getVelocity( xdot);
          
        	//compute Forces
        	computeForces(); 
	    	getForce(f);
	        	
	    	// Calculate rhs
	    	rhs.set(stepSize,f);// h*f
	    	
	    	for (Particle p: particles) {
	    		rhs.add(2*p.index,stepSize*p.df.x);
	    		rhs.add(2*p.index+1,stepSize*p.df.y);
	    	}
	    	
	    	// initialize deltaxdot
	    	deltaxdot.zero();
	      	
			//solve A deltaxdot = rhs
	    	cgMTJ.solve(this, rhs, deltaxdot, cgIterations.getValue()); 
	    	
	    	//update velocity
	    	xdot.add(deltaxdot);	
	    	xz.set(xdot);
           	setVelocity(xdot);
           	

        }
              
        // Finally, update the positions using the velocity at the next time step
        for ( Particle p : particles ) { 
        	
            if ( p.pinned ) continue;
            // symplectic Euler 
            tmp.scale( h, p.v );
            p.p.add( tmp );
            p.f.set(0,0);
          

        }
        
      processFracture();
                        
        time = time + h;
        return resultOK;
    }
    
    
    
	/**
	 * Comptues product with A = M - h (alpha M + beta K) - h^2 K
	 * But note there is an extra contribution for collision processing
	 */
	@Override
    public void mult(Vector v, Vector Av) {

		// TODO: Objective 2: implicit integration matrix multiply
		
    	// initialize 
    	int i = 0;
    	for (Particle  P : particles) {
    	//set the particle dx parameter
		P.dx.x = stepSize*v.get(i);
		P.dx.y = stepSize*v.get(i+1);
		i+=2;
    	//initialize the particle df to zero
    	P.df.set(0,0);
    	}
    	
    	// compute Kdx where dx=h*v  for each particle
     	for (FEMTriangle triangle : femSprings) 
     		triangle.computeForceDifferentials();
     	double c1 = useRAlpha.getValue() ? 1.0-RayleighAlpha.getValue()*stepSize :1.0;
     	double c2 = useRBeta.getValue() ? -RayleighBeta.getValue()-stepSize :-stepSize;
        for ( Particle p : particles ) {
        	
        	Av.set(2*p.index, c1*p.mass*v.get(2*p.index));
        	Av.set(2*p.index+1, c1*p.mass*v.get(2*p.index+1));
        	Av.add(2*p.index, c2*p.df.x);
        	Av.add(2*p.index+1, c2*p.df.y);
        }
       	
    	
    	
		for (Collision c : collisionList) {
			c.compute_dvDotResponse(Av, v, stepSize);
		}	
	}
	
	
	
	@Override
	public void filter(Vector v) {
		int i = 0;
		for ( Particle p : particles ) {
			if ( p.pinned ){
            	v.set( i+0, 0 );
            	v.set( i+1, 0 );
            }
			i += 2;
        }
	}

    /**
     * Computes all forces acting on the degrees of freedom in the system.
     */
    private void computeForces() {
    	// gravity
        for ( Particle p : particles ) {
            p.f.set( 0, useg.getValue() ? g.getValue() * p.mass : 0 );                      
        }
        
        for (Collision c : collisionList) {
        	c.applyForce( stepSize, implicit.getValue() );
        }
        
        if ( useMouseSpring ) {
        	mouseSpring.apply();
        }

        // TODO: Objective 1: apply the forces of the FEMSprings
       
        for (FEMTriangle triangle : femSprings) {
        	
        	triangle.applyForce();
    	
        }
        		
    	
        // TODO: Objective 2: explicit parts of the Rayleigh damping computation are needed here
        // use the computeForceDifferentials method of the the FEMTriangles, but note that you 
        // must set the particle dx parameter first, and initialize the particle df to zero as
        // differentials will be accumulated for all triangles adjacent to the particle.

	         for (Particle  P : particles) {
	        	//set the particle dx parameter
	    		P.dx.x =stepSize*P.v.x;
	    		P.dx.y = stepSize*P.v.y;
	        	//initialize the particle df to zero
	        	P.df.set(0,0);	
	        }
	        
	         // compute force differentials for each triangle
	    	for (FEMTriangle triangle : femSprings) 
	        	triangle.computeForceDifferentials();
	        	
	        
	    	
	    	// Add Rayleigh damping force  
	    	for (Particle p : particles) {
	    		
	    		p.f.x+=(useRBeta.getValue()? RayleighBeta.getValue()/stepSize:0.0)*p.df.x;
	    		p.f.y+=(useRBeta.getValue()? RayleighBeta.getValue()/stepSize:0.0)*p.df.y;
	    		
	    		p.f.x+=useRAlpha.getValue()? RayleighAlpha.getValue()*p.mass*p.v.x:0.0;
	    		p.f.y+=useRAlpha.getValue()? RayleighAlpha.getValue()*p.mass*p.v.y:0.0;

	    	}
    	  

    }
    
    /**
     * Process the fracture based on particles using the separation tensors.
     * The largest eigenvector that exceeds the toughness determines the 
     * line of separation. 
     */
    void processFracture() {
        // compute the separation tensors at each node
        for ( Particle p : particles ) {
        	if ( p.tris.size() == 0 ) continue;
        	p.computeSeparationTensor();
          
        }
        
        // TODO: Objective 4: process fracture
        // Note that you can use identifyBoundaries() as a slow way to update the border edges
        // after modifying the topology.

        double max = 0.0;
        int m_index=-1;
        int [] n_dgr = new int[particles.size()+1];
        for (Particle p : particles) {
        	n_dgr[p.index]=count_p_dgr(p.tris, p.index)-p.tris.size();
         	if (p.pinned ||  p.tris.size()<2)  continue;
           	max= Math.max(max,Math.max(p.separationTensor.ev1,p.separationTensor.ev2));
         	if (max == p.separationTensor.ev1 || max == p.separationTensor.ev2) 
         		m_index = p.index;
         }
        
        if ( m_index!=-1 && max > toughness.getValue() && particles.get(m_index).tris.size()>=1 ) {	
	        ArrayList<FEMTriangle> Clone_tris_1 = new ArrayList<FEMTriangle>(particles.get(m_index).tris );
	        int[] inds;
        	Vector2d normal = new Vector2d();
        	normal.set(max == particles.get(m_index).separationTensor.ev1 ? particles.get(m_index).separationTensor.v1 :  particles.get(m_index).separationTensor.v2 );
        	normal.normalize();
        	Particle Clone_p = new Particle(particles.get(m_index));
        	Clone_p.index=particles.size();
        	Clone_p.mass = particles.get(m_index).mass;

        	Split(m_index, normal, max, Clone_p,particles.get(m_index).tris);
        	ArrayList<FEMTriangle> Clone_tris = new ArrayList<FEMTriangle>(particles.get(m_index).tris );
        	
        	
	        if ( Clone_p.tris.size()!=0  ) 
			{  
	        	if(count_p_dgr(Clone_p.tris,particles.size())-Clone_p.tris.size()>=2){
	        		
	        	}

	        	else {
    	        particles.add(Clone_p);
	        	
    	        particles.get(m_index).tris.clear();
    	        particles.get(m_index).fplus.clear();
    	        particles.get(m_index).fminus.clear();
	        	
	    		for(FEMTriangle tri : Clone_tris) {
	    			
	    			if (tri.A.index==m_index) 
	    				
	    				tri.Ai= particles.get(m_index).addTriangle(tri);

	    			if (tri.B.index==m_index) 
    					tri.Bi=particles.get(m_index).addTriangle(tri);

    				if (tri.C.index==m_index) 
    					tri.Ci= particles.get(m_index).addTriangle(tri);
	    			
	    		}
	    		
	        	identifyBoundaries();
	        	}   		    		
			
			}
	        
        }
        
       } 
      
    
  
/**
 * 
 * @param m_index
 * @param max
 * @param n_dgr 
 * @return
 */
	private boolean Test_hinge_creation(ArrayList<FEMTriangle> tris,int m_index, double max, int[] n_dgr) {
		// TODO Auto-generated method stub
		Vector2d tmp = new Vector2d();
		
		Vector2d normal = new Vector2d();
		n_dgr[particles.size()]=0;
		
		for (FEMTriangle tri : tris) {
			normal.set(max == particles.get(m_index).separationTensor.ev1 ? particles.get(m_index).separationTensor.v1 :  particles.get(m_index).separationTensor.v2 );
        	normal.normalize();      	
        	tmp.set(particles.get(m_index).p);  
    		Point2d Center = new Point2d(tri.A.p.x+tri.B.p.x+tri.C.p.x,tri.A.p.y+tri.B.p.y+tri.C.p.y);
    		Center.scale(1.0/3.0);
    		tmp.sub(Center);
    		FEMTriangle T = null;
    		int i=0;
    		int n=n_dgr[m_index];
			T=adg_tri(tris,tri,m_index);
			if(T==null) return true;

    		if(normal.dot(tmp)>0) {	
    			i+=1;
    			tmp.set(particles.get(m_index).p);  
        		Center = new Point2d(T.A.p.x+T.B.p.x+T.C.p.x,T.A.p.y+T.B.p.y+T.C.p.y);
        		Center.scale(1.0/3.0);
        		tmp.sub(Center);		
    	    	if (normal.dot(tmp)>0) { 
    	    		n_dgr[particles.size()]=1;
    	    		
    	    	}
    	    	else {
    	    			
    	    			n_dgr[m_index]-=1;
    	    			n_dgr[particles.size()]+=1;
    	    		
    	    		
    	    	}
    	    }else {
    	    	
    	    }
    		
    		
		}
		if(n_dgr[m_index]>=2)return true;	
    	if(n_dgr[particles.size()]>=2)return true;	
	 return false;
	}

	private int get_ind(FEMTriangle t, FEMTriangle tri, int m_index) {
	// TODO Auto-generated method stub
		if(t.A.index == m_index) {
				if(t.B.index == tri.A.index || t.C.index == tri.A.index)return tri.A.index;
				if(t.B.index == tri.B.index || t.C.index == tri.B.index)return tri.B.index;
				if(t.B.index == tri.C.index || t.C.index == tri.C.index)return tri.C.index;

		}
		if(t.B.index == m_index) {
			if(t.A.index == tri.A.index || t.C.index == tri.A.index)return tri.A.index;
			if(t.A.index == tri.B.index || t.C.index == tri.B.index)return tri.B.index;
			if(t.A.index == tri.C.index || t.C.index == tri.C.index)return tri.C.index;

		}
		if(t.C.index == m_index) {
		if(t.B.index == tri.A.index || t.A.index == tri.A.index)return tri.A.index;
		if(t.B.index == tri.B.index || t.A.index == tri.B.index)return tri.B.index;
		if(t.B.index == tri.C.index || t.A.index == tri.C.index)return tri.C.index;

		}
		return -1;
		
}

	private FEMTriangle adg_tri(ArrayList<FEMTriangle> tris, FEMTriangle tri, int m_index) {
	// TODO Auto-generated method stub
		for (FEMTriangle T : tris) { 
    			if (T == tri) continue;
    			if(T.A.index == m_index) {
    				if(tri.A.index == m_index)
    					if(T.B.index == tri.C.index ||T.C.index == tri.C.index ||T.B.index == tri.B.index || T.C.index == tri.B.index)return T;
    				if(tri.B.index == m_index)
    					if(T.B.index == tri.C.index ||T.C.index == tri.C.index ||T.B.index == tri.A.index || T.C.index == tri.A.index)return T;
    				if(tri.C.index == m_index)
    					if(T.B.index == tri.A.index ||T.C.index == tri.A.index ||T.B.index == tri.B.index || T.C.index == tri.B.index)return T;
    			}
    			if(T.B.index == m_index) {
    				if(tri.A.index == m_index)
    					if(T.A.index == tri.C.index ||T.C.index == tri.C.index ||T.A.index == tri.B.index || T.C.index == tri.B.index)return T;
    				if(tri.B.index == m_index)
    					if(T.A.index == tri.C.index ||T.C.index == tri.C.index ||T.A.index == tri.A.index || T.C.index == tri.A.index)return T;
    				if(tri.C.index == m_index)
    					if(T.A.index == tri.A.index ||T.C.index == tri.A.index ||T.A.index == tri.B.index || T.C.index == tri.B.index)return T;
    			}
    			if(T.C.index == m_index) {
    				if(tri.A.index == m_index)
    					if(T.B.index == tri.C.index ||T.A.index == tri.C.index ||T.B.index == tri.B.index || T.A.index == tri.B.index)return T;
    				if(tri.B.index == m_index)
    					if(T.B.index == tri.C.index ||T.A.index == tri.C.index ||T.B.index == tri.A.index || T.A.index == tri.A.index)return T;
    				if(tri.C.index == m_index)
    					if(T.B.index == tri.A.index ||T.A.index == tri.A.index ||T.B.index == tri.B.index || T.A.index == tri.B.index)return T;
    			}

    		
		}
	return null;
}

	/**
     * Split the list of elements
     * @param tris
     * @param m_index
     * @param normal
     * @param max
     * @param clone_p 
     */
private void Split( int m_index, Vector2d normal, double max, Particle clone_p, ArrayList<FEMTriangle> tris) {
		// TODO Auto-generated method stub
	Vector2d tmp = new Vector2d(); 

    for (FEMTriangle Tri : particles.get(m_index).tris) {
    	tmp.set(particles.get(m_index).p);  
		Point2d Center = new Point2d(Tri.A.p.x+Tri.B.p.x+Tri.C.p.x,Tri.A.p.y+Tri.B.p.y+Tri.C.p.y);
		Center.scale(1.0/3.0);
		tmp.sub(Center); 
		
		if (normal.dot(tmp)>0.0) {
		    	if (Tri.A.index == m_index) {

		    		Tri.A=clone_p;
		    		Tri.Ai= clone_p.addTriangle(Tri);
					
		    	}
		    	if (Tri.B.index == m_index) { 

		    		Tri.B=clone_p;
		    		Tri.Bi= clone_p.addTriangle(Tri);
		    	
		    		}
		    	if (Tri.C.index == m_index) {

		    		Tri.C=clone_p;
		    		Tri.Ci= clone_p.addTriangle(Tri);
					
		    	}
			
		    }
		
		    else {	
		       	continue;
		    }
     	
	}
    
	
	
		
	}

/**
 * count the degree of a particle in the list of triangles containing the node
 * @param tris
 * @return
 */
    private int count_p_dgr(ArrayList<FEMTriangle> tris,int index) {
		// TODO Auto-generated method stub
    	int n=0;
    	for (int i= 0; i<tris.size(); i++ ) { 
    		for (int j=i+1 ;j<tris.size();j++ ) {		
    	
    			if (tris.get(j) == tris.get(i)) continue;
    					
    			if (tris.get(j).A.index==index) {
	    			if (tris.get(i).A.index==index) {
	    				if (tris.get(i).B== tris.get(j).B|| tris.get(i).C== tris.get(j).C||tris.get(i).B== tris.get(j).C|| tris.get(i).C== tris.get(j).B)
	    				n+=1;	
	    			}
	    			if (tris.get(i).B.index==index) {
	    				if (tris.get(i).A== tris.get(j).B|| tris.get(i).C== tris.get(j).C||tris.get(i).A== tris.get(j).C|| tris.get(i).C== tris.get(j).B)
		    				n+=1;
	    			}
	    			if (tris.get(i).C.index==index) {
	    				if (tris.get(i).B== tris.get(j).B|| tris.get(i).A== tris.get(j).C||tris.get(i).B== tris.get(j).C|| tris.get(i).A== tris.get(j).B)
		    				n+=1;
	    			}
    			}
    			else {
    			if (tris.get(j).B.index==index) {
	    			if (tris.get(i).A.index==index) {
	    				if (tris.get(i).B== tris.get(j).A|| tris.get(i).C== tris.get(j).C||tris.get(i).B== tris.get(j).C|| tris.get(i).C == tris.get(j).A)
	    				n+=1;	
	    			}
	    			if (tris.get(i).B.index==index) {
	    				if (tris.get(i).A== tris.get(j).A|| tris.get(i).C== tris.get(j).C||tris.get(i).A== tris.get(j).C|| tris.get(i).C == tris.get(j).A)
		    				n+=1;	
	    			}
	    			if (tris.get(i).C.index==index) {
	    				if (tris.get(i).B== tris.get(j).A|| tris.get(i).A== tris.get(j).C||tris.get(i).B== tris.get(j).C|| tris.get(i).A == tris.get(j).A)
		    				n+=1;	
	    			}
    			}
    			}
    			if (tris.get(j).C.index==index) {
	    			if (tris.get(i).A.index==index) {
	    				if (tris.get(i).B== tris.get(j).A|| tris.get(i).C== tris.get(j).B||tris.get(i).B== tris.get(j).B|| tris.get(i).C == tris.get(j).A)
	    				n+=1;		
	    			}
	    			if (tris.get(i).B.index==index) {
	    				if (tris.get(i).A== tris.get(j).A|| tris.get(i).C== tris.get(j).B||tris.get(i).A== tris.get(j).B|| tris.get(i).C == tris.get(j).A)
		    				n+=1;
	    			}
	    			if (tris.get(i).C.index==index) {
	    				if (tris.get(i).B == tris.get(j).A|| tris.get(i).A== tris.get(j).B||tris.get(i).B== tris.get(j).B|| tris.get(i).A == tris.get(j).A)
		    				n+=1;
	    			}
    			}
    			
    		}

    	}
    	n+=(tris.size()-n)*2;
		return n;
	}

    
    
    
	DenseVector xdot;
    DenseVector xz;
    DenseVector deltax;
    DenseVector deltaxdot;
    DenseVector f;
    DenseVector rhs;
    ConjugateGradientMTJ cgMTJ;

    /**
     * Initializes variables for backward Euler integration.
     * You may not actually need all these.  Rename them as you like
     */
    public void init() {
        int N = particles.size();
        if ( xdot == null || xdot.size() != 2*N ) {
        	xdot = new DenseVector( 2*N );
        	xz = new DenseVector( 2*N );
        	deltax = new DenseVector( 2*N );
        	deltaxdot = new DenseVector( 2*N );
        	f = new DenseVector( 2*N );
            rhs = new DenseVector( 2*N );
        	cgMTJ = new ConjugateGradientMTJ( 2*N );
            cgMTJ.setFilter(this);
            
        }      
   
    }
    
    
    /**
     * Fills in the provided vector with the particle velocities.
     * @param xd
     */
    private void getVelocity( DenseVector xd ) {
    	int j = 0;
        for ( Particle p : particles ) {
            if( p.pinned ) {
                xd.set( j+0, 0 );
                xd.set( j+1, 0 );
            } else {
                xd.set( j+0, p.v.x );
                xd.set( j+1, p.v.y );
            }
            j += 2;
        }       
    }

    private void setVelocity( DenseVector xdot ) {
    	int j = 0;
        for ( Particle p : particles ) {
            if( p.pinned ) {
                p.v.set(0,0);
            } else {
                p.v.x = xdot.get(j+0);
                p.v.y = xdot.get(j+1);
            }
            j += 2;
        }
    }
    
    private void getForce( DenseVector f ) {
    	int j = 0;
        for ( Particle p : particles ) {
        	f.set( j+0, p.f.x );
        	f.set( j+1, p.f.y );
        	j += 2;
        }
    }
    
  
    /**
     * set Positions 
     * @param x
     */
    private void setPosition ( DenseVector x ) {
    	int j = 0;
        for ( Particle p : particles ) {
    
                p.p.x = x.get(j+0);
                p.p.y = x.get(j+1); 
            j += 2;
        }
    }
    
    
    public void applyTempVelocities() {
    	for(Particle p : particles)	{
    		p.v.set(p.vTmp);
    	}
    }
    
    public void updateTempVelocities() {
    	for(Particle p : particles) {
    		p.vTmp.set(p.v);
    	}
    }
    /**
     * Creates a new particle and adds it to the system
     * @param x
     * @param y
     * @param vx
     * @param vy
     * @return the new particle
     */
    public Particle createParticle( double x, double y, double vx, double vy ) {
        Particle p = new Particle( x, y, vx, vy );
        particles.add( p );
        return p;
    }
            
    @Override
    public void init(GLAutoDrawable drawable) {
        // do nothing
    }

    @Override
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();

        // draw the interaction spring
        if ( useMouseSpring ) {
        	mouseSpring.display( drawable );        	
        }

        // draw the boundaries of the mesh (i.e., used in collision detection
        if ( drawBoundaryEdges.getValue() ) {
	        gl.glColor4d( 0,0,0, 1 );
	        gl.glLineWidth( 2 );
	        gl.glBegin( GL.GL_LINES );
	        for ( Edge e : collidableEdges ) {
	        	gl.glVertex2d( e.p1.p.x, e.p1.p.y );
	        	gl.glVertex2d( e.p2.p.x, e.p2.p.y );
	        }
	        gl.glEnd();
        }
        
        // draw the triangles, and stress
        gl.glLineWidth( 1 );
        float alpha = transparency.getFloatValue();
        
        if ( drawElementBoundaries.getValue() ) {
        	for ( FEMTriangle f : femSprings ) {
        		f.displayElementBoundaries(drawable, alpha);
        	}
        }
        
        double s1 = strainEigVecScale.getValue();
        double s2 = stressEigVecScale.getValue();
        for ( FEMTriangle f : femSprings ) {
           f.display(drawable, alpha);
           if ( drawStrainTensor.getValue() ) {
        	   f.displayStrain( drawable, s1 );
           }
           if ( drawStressTensor.getValue() ) {
            	f.displayStress( drawable, s2 );
            }
        }

        // draw the separation tensors
        double sts = stScale.getValue();        
        if ( drawSeparationTensor.getValue() ) {
        	for ( Particle p : particles ) {
        		p.drawSeparationTensor(drawable, sts);
        	}
        }
        
        // draw the particles
        if ( drawParticles.getValue() ) {
            gl.glPointSize( pointSize.getFloatValue() );
            gl.glBegin( GL.GL_POINTS );
            for ( Particle p : particles ) {
                // transparency is used to get smooth edges on the particles
                alpha = 1;//  0.75;
                if ( p.pinned ) {
                    gl.glColor4d( 1, 0, 0, alpha );
                } else {
                	gl.glColor4d( 0,0,0, alpha );//gl.glColor4d( 0, 0.95,0, alpha );
                }
                gl.glVertex2d( p.p.x, p.p.y );
            }
            gl.glEnd();
        }
        
        if ( drawCollisions.getValue() ) {
        	for ( Collision c : collisionList ) {
        		c.display( drawable, 1.0 );
        	}
        }
        if ( drawCollisionBondary.getValue() ) {
        	for ( Collision c : collisionList ) {
        		c.drawCollisionBoundary(drawable);
        	}
        }
        
        
    }    
    
    // visualizaiton related parameters
    
    BooleanParameter drawParticles = new BooleanParameter( "draw particles", true ) ;
    BooleanParameter drawCollisions = new BooleanParameter( "draw collisions", true ) ;
    BooleanParameter drawCollisionBondary = new BooleanParameter( "draw collision bondary", false ) ;
    DoubleParameter pointSize = new DoubleParameter("point size", 2, 1, 25);
    DoubleParameter stressEigVecScale = new DoubleParameter("stress eigenvector viz scale", 0.1, 1e-3, 1e3 );
    DoubleParameter strainEigVecScale = new DoubleParameter("strain eigenvector viz scale", 1e3, 1, 1e6 );
    DoubleParameter stScale = new DoubleParameter("separation tessor viz scale", 1e-7, 1e-10, 1 );
    BooleanParameter drawSeparationTensor = new BooleanParameter("draw separation tensor", false );
    BooleanParameter drawStrainTensor = new BooleanParameter("draw strain tensor", false );
    BooleanParameter drawStressTensor = new BooleanParameter("draw stress tensor", false );
    BooleanParameter drawElementBoundaries = new BooleanParameter( "draw element boundaries", false );
    BooleanParameter drawBoundaryEdges = new BooleanParameter( "draw boundary edges", true );
    DoubleParameter transparency = new DoubleParameter("triangle transparency", 0.5, 0, 1 );
    JTextArea comments = new JTextArea("<comments>");
    BooleanParameter showCommentsAndParameters = new BooleanParameter("show comments and parameters", true );
    
    // simulation related parameters
    
    /** 
     * Irving presents a solution for inverting elements, thus this makes an interesting
     * alternative material (not part of this assignment)
     */
    BooleanParameter useIrving2004 = new BooleanParameter( "use Irving 2004 (otherwise StVK)", false );
    
    BooleanParameter implicit = new BooleanParameter( "Implicit integration", true );
    IntParameter newtonIterations = new IntParameter( "Newton root solve iterations", 1, 1, 20 );
    IntParameter cgIterations = new IntParameter( "CG solve iterations", 100, 20, 300 );
    BooleanParameter doCollisionDetection = new BooleanParameter( "do collisions", true ) ;
    BooleanParameter useg = new BooleanParameter( "use gravity", true );
    DoubleParameter g = new DoubleParameter( "gravity", 98, 0.01, 1000 );
    DoubleParameter mouseSpringk = new DoubleParameter( "mouse spring k", 1000, 1, 1e5 );
    DoubleParameter PoissonRatio  = new DoubleParameter( "Poisson Ratio", 0.3, 0, .4999 );
    DoubleParameter YoungModulus = new DoubleParameter( "YoungModulus", 50000, 1, 1e10 );
    DoubleParameter toughness = new DoubleParameter("material toughness", 1e5, 1, 1e8 );
    DoubleParameter RayleighAlpha = new DoubleParameter("Raleigh alpha", 1e-1, 1e-3, 1e3 );
    DoubleParameter RayleighBeta = new DoubleParameter("Raleigh beta", 1e-2, 1e-3, 1e3 );
    BooleanParameter useRAlpha = new BooleanParameter("use Raleigh alpha", true );
    BooleanParameter useRBeta = new BooleanParameter("use Raleigh beta", true );
    
    /** coefficient for viscous resolution of collision */
	DoubleParameter collisionViscousCoefficient = new DoubleParameter( "collision viscous coefficient", 5e7, 1, 5e8 );
    /** relative normal velocity threshold for using damping based collision */
	DoubleParameter relNormalVelThresh = new DoubleParameter( "relative normal velocity threshold", 100, 1e-3, 1e3 );
	/** Could be area weighted... was (2.0 * area + 50); */
	DoubleParameter collisionRepulsion = new DoubleParameter( "collision repulsion (buggy)", 0, 1e-3, 1e3 ); 
    
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        
        VerticalFlowPanel vfp0 = new VerticalFlowPanel();
        vfp0.setBorder( new TitledBorder("Viewing Parameters" ) );
        
        HorizontalFlowPanel hfp0 = new HorizontalFlowPanel();
        hfp0.add( drawParticles.getControls() );
        hfp0.add( pointSize.getSliderControls(false) );
        vfp0.add( hfp0.getPanel() );
        
        vfp0.add( drawCollisions.getControls() );
        vfp0.add( drawCollisionBondary.getControls() );
        
        HorizontalFlowPanel hfp1 = new HorizontalFlowPanel();
        hfp1.add( drawStrainTensor.getControls() );
        hfp1.add( strainEigVecScale.getSliderControls(true));
        vfp0.add( hfp1.getPanel() );
        
        HorizontalFlowPanel hfp2 = new HorizontalFlowPanel();
        hfp2.add( drawStressTensor.getControls() );
        hfp2.add( stressEigVecScale.getSliderControls(true));
        vfp0.add( hfp2.getPanel() );
        
        HorizontalFlowPanel hfp3 = new HorizontalFlowPanel();
        hfp3.add( drawSeparationTensor.getControls() );
        hfp3.add( stScale.getSliderControls(true));
        vfp0.add( hfp3.getPanel() );
        
        vfp0.add( drawElementBoundaries.getControls() );
        vfp0.add( drawBoundaryEdges.getControls() );
        vfp0.add( transparency.getSliderControls(false) );
        vfp0.add( comments );
        vfp0.add( showCommentsAndParameters.getControls() );
        CollapsiblePanel cp0 = new CollapsiblePanel( vfp0.getPanel() );
        cp0.collapse();
        vfp.add( cp0 );
                
        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder( new TitledBorder("Simulation Parameters"));
        vfp1.add( doCollisionDetection.getControls() );
        vfp1.add( implicit.getControls() );
        vfp1.add( newtonIterations.getControls() );
        vfp1.add( cgIterations.getControls() );
        
        HorizontalFlowPanel hfpg = new HorizontalFlowPanel();
        hfpg.add( useg.getControls() );
        hfpg.add( g.getSliderControls(true) );
        vfp1.add( hfpg.getPanel() );

        vfp1.add( YoungModulus.getSliderControls(true));
        vfp1.add( PoissonRatio.getSliderControls(false));
        vfp1.add( mouseSpringk.getSliderControls(true));
        vfp1.add( toughness.getSliderControls(true));
        
        HorizontalFlowPanel hfp4 = new HorizontalFlowPanel();
        hfp4.add( useRAlpha.getControls() );
        hfp4.add( RayleighAlpha.getSliderControls(true));
        vfp1.add( hfp4.getPanel() );
        
        HorizontalFlowPanel hfp5 = new HorizontalFlowPanel();
        hfp5.add( useRBeta.getControls() );
        hfp5.add( RayleighBeta.getSliderControls(true));
        vfp1.add( hfp5.getPanel() );
        
        vfp1.add( useIrving2004.getControls() );
        vfp1.add( collisionViscousCoefficient.getSliderControls(true) );
        vfp1.add( collisionRepulsion.getSliderControls(true) );
        vfp1.add( relNormalVelThresh.getSliderControls(true) );
        
        CollapsiblePanel cp1 = new CollapsiblePanel( vfp1.getPanel() );
        cp1.collapse();
        vfp.add( cp1 );
        
        return vfp.getPanel();        
    }
    
    @Override
    public String toString() {
    	DecimalFormat df = new DecimalFormat("0.000");
        String s = "particles = " + particles.size() + "\n" + "time = " + df.format(time) + "\n";
        s += comments.getText() + "\n";                
        return s;
    }
    
}
