package fractureW20;

import java.util.ArrayList;
import java.util.Comparator;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector2d;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

public class Collision {
	
	/** Triangle 1 in collision */
	FEMTriangle t1;
	
	/** Triangle 2 in collision */
	FEMTriangle t2;

	/** Area of the overlap */
	double area;

	/** Centroid of the overlap between two triangles */
	Point2d centroid = new Point2d();
	
	/** Direction of normal force acting between two triangles, unit length */
	Vector2d normal = new Vector2d(); 	

	/** coefficient controlling the magnitude of viscous collision */
	double coefficent_viscous;
	/** coefficient controlling the magnitude of repulsive collision */
	double coefficent_repulsive;
	/** low relative velocity threshold */
	double relNormalVelThres;
	/** This is used for testing for zeros to avoid numerical errors. */
	double eps;
		
	/**
	 * List of points on collision boundary, exclusively used for debugging intersections
	 * (this could be removed to speed things up slightly)
	 */
	ArrayList<Point2d> collisionBoundary = new ArrayList<Point2d>();
	
	/** viscous response */
	DenseMatrix response_viscous = new DenseMatrix(12, 12);
	
	/** spring response */
	DenseMatrix response_repulse = new DenseMatrix(12, 12);
	
	DenseVector velocities = new DenseVector(12);
	
	DenseVector positions = new DenseVector(12);
	
	DenseVector forces_viscous = new DenseVector(12);
	
	DenseVector forces_repulsive = new DenseVector(12);
	
	/**
	 * Creates a new collision object
	 * @param t1
	 * @param t2
	 * @param centroid
	 * @param normalT1
	 * @param normalT2
	 * @param collisionBoundary
	 * @param eps
	 */
	public Collision(FEMTriangle t1, FEMTriangle t2, Point2d centroid, double area, Vector2d normal, ArrayList<Point2d> collisionBoundary, double eps) {
		this.t1 = t1;
		this.t2 = t2;
		this.centroid.set( centroid );
		this.collisionBoundary = collisionBoundary;
		this.eps = eps;
		this.area = area; 
		this.normal.set( normal );
		
		response_viscous.zero();
		response_repulse.zero();
		
		// Default collision respose parameters will typically be adjusted just
		// after this constructor is called.
		coefficent_viscous = 50000000;
		relNormalVelThres = 100;
		coefficent_repulsive = (2.0 * area + 50); // negative or positive?
				
		Vector2d centroidT1 = new Vector2d(t1.A.p);
		centroidT1.add(t1.B.p);
		centroidT1.add(t1.C.p);
		Vector2d centroidT2 = new Vector2d(t2.A.p);
		centroidT2.add(t2.B.p);
		centroidT2.add(t2.C.p);
		centroidT2.sub(centroidT1);
				
		if (centroidT2.length() > eps) {
			coefficent_repulsive /= centroidT2.length();
		} else {
			coefficent_repulsive = 0;
		}
		
		coefficent_repulsive *= -1;
		
		int index = 0;
		velocities.set(index++, t1.A.v.x);
		velocities.set(index++, t1.A.v.y);
		velocities.set(index++, t1.B.v.x);
		velocities.set(index++, t1.B.v.y);
		velocities.set(index++, t1.C.v.x);
		velocities.set(index++, t1.C.v.y);
		velocities.set(index++, t2.A.v.x);
		velocities.set(index++, t2.A.v.y);
		velocities.set(index++, t2.B.v.x);
		velocities.set(index++, t2.B.v.y);
		velocities.set(index++, t2.C.v.x);
		velocities.set(index++, t2.C.v.y);
		
		index = 0;
		positions.set(index++, t1.A.p.x);
		positions.set(index++, t1.A.p.y);
		positions.set(index++, t1.B.p.x);
		positions.set(index++, t1.B.p.y);
		positions.set(index++, t1.C.p.x);
		positions.set(index++, t1.C.p.y);
		positions.set(index++, t2.A.p.x);
		positions.set(index++, t2.A.p.y);
		positions.set(index++, t2.B.p.x);
		positions.set(index++, t2.B.p.y);
		positions.set(index++, t2.C.p.x);
		positions.set(index++, t2.C.p.y);
		
	}
	
	/**
	 * Draws collision points and normals
	 * @param drawable
	 * @param scale how much to scale the normal vector
	 */
	public void display( GLAutoDrawable drawable, double scale ) {
		GL2 gl = drawable.getGL().getGL2();
		gl.glColor3f(1, 0, 0);
		gl.glPointSize(4f);
		gl.glBegin( GL.GL_POINTS );
		gl.glVertex2d( centroid.x, centroid.y );
		gl.glEnd();
		gl.glBegin( GL.GL_LINES );
		gl.glColor3f(0.95f, 0.3f, 0.95f);
		gl.glVertex2d( centroid.x, centroid.y );
		gl.glVertex2d( centroid.x + normal.x * scale, centroid.y + normal.y * scale );
		gl.glEnd();
	}
	
	/**
	 * Draws additional collision detection information
	 * such as the whole boundry and forces on the nodes
	 * @param drawable
	 */
	public void drawCollisionBoundary( GLAutoDrawable drawable ) {
		// For visualization only
		DenseVector forceStart = new DenseVector(12);
		DenseVector forceEnd = new DenseVector(12);
		
		double c1 = 0; // show viscous
		double c2 = 10; // show repulsive
						
		int index = 0;
		forceStart.set(index++, t1.A.p.x);
		forceStart.set(index++, t1.A.p.y);
		forceStart.set(index++, t1.B.p.x);
		forceStart.set(index++, t1.B.p.y);
		forceStart.set(index++, t1.C.p.x);
		forceStart.set(index++, t1.C.p.y);
		forceStart.set(index++, t2.A.p.x);
		forceStart.set(index++, t2.A.p.y);
		forceStart.set(index++, t2.B.p.x);
		forceStart.set(index++, t2.B.p.y);
		forceStart.set(index++, t2.C.p.x);
		forceStart.set(index++, t2.C.p.y);
		
		forceEnd.set(1, forceStart);
		forceEnd.add(c1 * 0.001, forces_viscous);
		forceEnd.add(c2 * 0.001, forces_repulsive);
		
		GL2 gl = drawable.getGL().getGL2();
			
		gl.glBegin( GL.GL_POINTS );
		gl.glColor3f(0, 0.9f, 0.95f);
		for (Point2d pp : collisionBoundary) {
			gl.glVertex2d( pp.x, pp.y );
		}
		gl.glEnd();
		
		gl.glBegin( GL.GL_LINE_LOOP );
		gl.glColor3f(0, 0.9f, 0.95f);
		for (Point2d pp : collisionBoundary) {
			gl.glVertex2d( pp.x, pp.y );
		}
		gl.glEnd();
		
		gl.glBegin( GL.GL_LINES );
		for (int i = 0; i < 12; i += 2) {
			float r = (i + 6) % 6 == 0 ? 1 : 0;
			float g = (i + 4) % 6 == 0 ? 1 : 0;
			float b = (i + 2) % 6 == 0 ? 1 : 0;
			gl.glColor3f(r, g, b);
			gl.glVertex2d( forceStart.get(i), forceStart.get(i + 1) );
			gl.glVertex2d( forceEnd.get(i), forceEnd.get(i + 1) );
		}
		gl.glEnd();
	}
		
	
	/**
	 * Computes the intersection between two line segments.
	 * @param p1 start first line segment
	 * @param p2 end first line segment
	 * @param q1 start second line segment
	 * @param q2 end second line segment
	 * @param eps threshold for dealing with tricky cases
	 * @param intersection will contain the intersection point if it exists, otherwise no change
	 * @return true on intersection, in which case intersection will be updated
	 */
	private static boolean lineLineIntersect( Point2d p1, Point2d p2, Point2d q1, Point2d q2, double eps, Point2d intersection ) {	
		Vector2d v1 = new Vector2d(p2);
		v1.sub(p1);
		Vector2d v2 = new Vector2d(q2);
		v2.sub(q1);
		
		double det = v2.x * v1.y - v1.x * v2.y;
		// Assume no intersection when lines are parallel.
		if ( Math.abs(det) < eps ) return false;
					
		double s, t;
		
		Vector2d pq = new Vector2d(p1);
		pq.sub(q1);
		
		s = (v1.y * pq.x - v1.x * pq.y) / det;
		
		if (s < 0 || s > 1) return false;
		
		if ( Math.abs(v1.x) >  Math.abs(v1.y) ) {
			t = (v2.x * s - pq.x) / v1.x;
		} else {
			t = (v2.y * s - pq.y) / v1.y;
		}
		
		if (t < 0 || t > 1) return false;
				
		if (t > s) {
			intersection.set( p1.x + v1.x * t, p1.y + v1.y * t );
			return true;
		}
		intersection.set( q1.x + v2.x * s, q1.y + v2.y * s );
		return true;
	}
	
	final static class MyComparator implements Comparator<Point2d> {
		Point2d centroid;
		Vector2d reference;
		public MyComparator( Point2d centroid, Vector2d reference ) {
			this.centroid = centroid;
			this.reference = reference;			
		}
	   public int compare(Point2d p1, Point2d p2) {
		   Vector2d v1 = new Vector2d(p1);
		   v1.sub(centroid);
		   v1.normalize();
		   
		   Vector2d v2 = new Vector2d(p2);
		   v2.sub(centroid);
		   v2.normalize();
		   
		   double theta1 = Math.atan2(reference.x * v1.y - reference.y * v1.x, reference.dot(v1));
		   theta1 = theta1 < 0 ? 2 * Math.PI + theta1 : theta1;
		   double theta2 = Math.atan2(reference.x * v2.y - reference.y * v2.x, reference.dot(v2));
		   theta2 = theta2 < 0 ? 2 * Math.PI + theta2 : theta2;
		   
		   return theta1 > theta2 ? 1 : -1;
	   }
	}
	
	/**
	 * Computes the centroid, and sorts the colList aroudn the centroid
	 * @param colList
	 * @param centroid
	 */
	private static void sortCClockwise( ArrayList<Point2d> colList, Point2d centroid ) {
		centroid.x = 0;
		centroid.y = 0;		
		for (Point2d p : colList) {
			centroid.add(p);
		}		
		centroid.scale(1.0 / colList.size());
		Vector2d reference = new Vector2d(colList.get(0));
		reference.sub(centroid);
		reference.normalize();		
		colList.sort( new MyComparator( centroid, reference ) );			   
	}
	
	/**
	 * @param list
	 * @param p
	 * @param eps
	 * @return true if any point in list is within eps of p
	 */
	private static boolean contains( ArrayList<Point2d> list, Point2d p, double eps ) {
		double avgDist = 0;
		for ( Point2d i : list ) {
			avgDist += i.distance(p);
		}
		avgDist /= list.size();
		for ( Point2d i : list ) {
			if ( i.distance(p) < eps * avgDist) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Performs collision detection 
	 * Creates and returns a Collision object if there is overlab between the two triangles
	 * @param t1
	 * @param t2
	 * @param eps
	 * @return null if no collision
	 */
	public static Collision collisionDetect( FEMTriangle t1, FEMTriangle t2, double eps ) {
		final ArrayList<Point2d> collisionList = new ArrayList<Point2d>();
		collisionList.clear();
		final Point2d centroid = new Point2d();
		final Vector2d normalT1 = new Vector2d();
		final Vector2d normalT2 = new Vector2d();
		final Vector2d normal = new Vector2d();
		
		// Check if the vertex of a triangle is inside another.
		if ( t1.isInside( t2.A.p, eps, null ) ) collisionList.add( t2.A.p );
		if ( t1.isInside( t2.B.p, eps, null ) ) collisionList.add( t2.B.p );
		if ( t1.isInside( t2.C.p, eps, null ) ) collisionList.add( t2.C.p );
		if ( t2.isInside( t1.A.p, eps, null ) ) collisionList.add( t1.A.p );
		if ( t2.isInside( t1.B.p, eps, null ) ) collisionList.add( t1.B.p );
		if ( t2.isInside( t1.C.p, eps, null ) ) collisionList.add( t1.C.p );
				
		// Check for intersection of edges
		final Point2d intersection = new Point2d();
		if ( lineLineIntersect(t1.A.p, t1.B.p, t2.A.p, t2.B.p, eps, intersection) ) collisionList.add( new Point2d( intersection ) );
		if ( lineLineIntersect(t1.A.p, t1.B.p, t2.B.p, t2.C.p, eps, intersection) ) collisionList.add( new Point2d( intersection ) );
		if ( lineLineIntersect(t1.A.p, t1.B.p, t2.C.p, t2.A.p, eps, intersection) ) collisionList.add( new Point2d( intersection ) );
		if ( lineLineIntersect(t1.B.p, t1.C.p, t2.A.p, t2.B.p, eps, intersection) ) collisionList.add( new Point2d( intersection ) );
		if ( lineLineIntersect(t1.B.p, t1.C.p, t2.B.p, t2.C.p, eps, intersection) ) collisionList.add( new Point2d( intersection ) );
		if ( lineLineIntersect(t1.B.p, t1.C.p, t2.C.p, t2.A.p, eps, intersection) ) collisionList.add( new Point2d( intersection ) );
		if ( lineLineIntersect(t1.C.p, t1.A.p, t2.A.p, t2.B.p, eps, intersection) ) collisionList.add( new Point2d( intersection ) );
		if ( lineLineIntersect(t1.C.p, t1.A.p, t2.B.p, t2.C.p, eps, intersection) ) collisionList.add( new Point2d( intersection ) );
		if ( lineLineIntersect(t1.C.p, t1.A.p, t2.C.p, t2.A.p, eps, intersection) ) collisionList.add( new Point2d( intersection ) );
				
		// Remove duplicates... a bit slow to do things this way.
		final ArrayList<Point2d> tempList = new ArrayList<Point2d>();
		for (Point2d c : collisionList) {
			if ( ! contains(tempList, c, eps) ) {
				tempList.add(c);
			}
		}		
	    collisionList.clear(); 
	    collisionList.addAll(tempList); 
	    	    
		// If number of intersection detected is less than 3, assume no collision.
		if (collisionList.size() < 3) {
			collisionList.clear();
			return null;
		}
				
		// Compute the centroid and sort the vertcies defining collision boundary in CCW order
		sortCClockwise( collisionList, centroid );
		
		Vector2d edge = new Vector2d(0, 0);
		double area = 0;
		normalT1.set(0, 0);
		normalT2.set(0, 0);
		
		centroid.x = 0;
		centroid.y = 0;
		double weightCentroid = 0;
		Point2d a = collisionList.get(0);
		for ( int i = 1; i <= collisionList.size(); i++ ) {
			Point2d b = collisionList.get( i % collisionList.size() );
			area += a.x * b.y - b.x * a.y;
			
			edge.sub(a, b);
			centroid.x += (a.x + b.x) * edge.length();
			centroid.y += (a.y + b.y) * edge.length();
			weightCentroid += edge.length();
			
			boolean edgeBelongsToT1 = t1.isSegmentOnTraingle(a, b, eps);
			boolean edgeBelongsToT2 = t2.isSegmentOnTraingle(a, b, eps);
			if ((edgeBelongsToT1 && edgeBelongsToT2) || !(edgeBelongsToT1 || edgeBelongsToT2)) {
				normalT1.x += edge.y;
				normalT1.y -= edge.x;
			}
			else if (edgeBelongsToT1) {
				normalT1.x += edge.y;
				normalT1.y -= edge.x;
			}
			else if (edgeBelongsToT2) {
				normalT2.x += edge.y;
				normalT2.y -= edge.x;
			}
			
			a = b;
		}
		area *= 0.5;
		area = Math.abs(area);
				
		centroid.scale(0.5 / weightCentroid);
				
		normal.add(normalT1, normalT2); // Should be zero.
		// if the intersection area is too small, or if there is no good normal information
		// then give up on this contact!
		if ( Math.abs(area) < eps * 10  || normal.length() > eps || (normalT1.length() + normalT2.length()) < 2 * eps) {
			return null;
		}
			
		normal.sub( normalT1, normalT2 );
		normal.normalize();
		
		//if (Double.isNaN(normal.x) || Double.isNaN(normal.y) || Double.isInfinite(normal.x) || Double.isInfinite(normal.y))
			//System.out.println("Problem");
		
		return new Collision( t1, t2, centroid, area, normal, collisionList, eps );
	}
	
	/**
	 * Computes the matrices needed for explicit and implicit damping based collision response
	 */
	public void updateCollisionResponse() {
		// Get the barycentric coordinates of centroid wrt triangle 1 and 2.
		final Point3d coords1 = new Point3d();
		final Point3d coords2 = new Point3d();
		boolean b1 = t1.isInside( centroid, eps, coords1 );
		boolean b2 = t2.isInside( centroid, eps, coords2 );
		assert( b1 && b2 );
		
		// Matrix of barycentric weights.
		DenseMatrix B = new DenseMatrix(2, 12);
		B.set(0, 0,  coords1.x); B.set(0, 1, 0); B.set(0, 2,  coords1.y); B.set(0, 3, 0); B.set(0,  4,  coords1.z ); B.set(0,  5, 0);
		B.set(0, 6, -coords2.x); B.set(0, 7, 0); B.set(0, 8, -coords2.y); B.set(0, 9, 0); B.set(0, 10, -coords2.z ); B.set(0, 11, 0);
		
		// Matrix of centroid, for repulsive response
		DenseMatrix C = new DenseMatrix(2, 12);
		C.set(0, 0, 1); C.set(0, 1, 0); C.set(0, 2, 1); C.set(0, 3, 0); C.set(0, 4, 1); C.set(0, 5, 0);
		C.set(0, 6, -1); C.set(0, 7, 0); C.set(0, 8, -1); C.set(0, 9, 0); C.set(0, 10, -1); C.set(0, 11, 0);
		
		for ( int i = 0; i < 12; i++) {
			if (i%2 == 0) {
				B.set(1, i, 0);
				C.set(1, i, 0);
			} else {
				B.set(1, i, B.get(0, i - 1));
				C.set(1, i, C.get(0, i - 1));
			}
		}
					
		// Negative transpose of barycentric weights
		DenseMatrix BT = new DenseMatrix(12, 2);
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 2; j++) {
				double element = B.get(j, i);
				if (Math.abs(element) > eps) {
					element = -element;
				}
				BT.set(i, j, element);
			}
		}
				
		DenseVector relativeVelocity = new DenseVector(2);
		B.mult(velocities, relativeVelocity);
			
		double nDotRelVel = relativeVelocity.get(0) * normal.x + relativeVelocity.get(1) * normal.y;
				
		DenseMatrix n = new DenseMatrix(2, 1);
		n.set(0, 0, normal.x);
		n.set(1, 0, normal.y);
		DenseMatrix nt = new DenseMatrix(1, 2);
		nt.set(0, 0, normal.x);
		nt.set(0, 1, normal.y);
		
		DenseMatrix temp1 = new DenseMatrix(12, 1);
		BT.mult(n, temp1);
		DenseMatrix temp2 = new DenseMatrix(12, 2);
		temp1.mult(nt, temp2);
		temp2.mult(C, response_repulse);
		temp2.mult(B, response_viscous);
		
		response_viscous.scale(coefficent_viscous);
		response_repulse.scale(coefficent_repulsive);
		
		// Unilateral response
		// We scale back the response if relative velocity is low and triangles are seprating.
		// If relative velocity is low
	 	if ( Math.abs(nDotRelVel) < relNormalVelThres ) {
			Vector2d centroidT1 = new Vector2d(t1.A.p);
			centroidT1.add(t1.B.p); centroidT1.add(t1.C.p);
			centroidT1.scale(1.0/3);
			centroidT1.sub(centroid); // Vector between centroid of triangle t1 and centroid of collision polygon.
						
			double nDotCentroid = centroidT1.dot(normal);
											
			// Two conditions for seperating triangle
			if (nDotCentroid > 0 && nDotRelVel > 0)
					response_viscous.scale(0);
			if (nDotCentroid <= 0 && nDotRelVel <= 0)
					response_viscous.scale(0);
					
		}
	}
	
	static private DenseVector dvLocal = new DenseVector(12);
	static private DenseVector dvDotResponseLocal = new DenseVector(12);
	static private int indexBuffer[] = new int[6];

	/**
	 * Computes the product of a given vector with the damping matrix for this collision.
	 * Av = Av + A * dv where A is the damping matrix.
	 * @param Av Output vector, to accumulate the result.
	 * @param dv the vector to be multiplied
	 * @param stepSize
	 */
	public void compute_dvDotResponse(Vector Av, Vector dv, double stepSize) {
		indexBuffer[0] = t1.A.index;
		indexBuffer[1] = t1.B.index;
		indexBuffer[2] = t1.C.index;
		indexBuffer[3] = t2.A.index;
		indexBuffer[4] = t2.B.index;
		indexBuffer[5] = t2.C.index;
		
		// Use a local view of only the combined 12 DOFs of each triangle.
		dvLocal.zero();
		dvDotResponseLocal.zero();
		
		for (int i = 0; i < 6; i++) {
			dvLocal.set(2 * i, dv.get(indexBuffer[i] * 2));
			dvLocal.set(2 * i + 1, dv.get(indexBuffer[i] * 2 + 1));
		}
		
		response_viscous.mult(dvLocal, dvDotResponseLocal);
		dvDotResponseLocal.scale(-stepSize);
				
		for (int i = 0; i < 6; i++) {
			Av.add(indexBuffer[i] * 2, dvDotResponseLocal.get(2 * i));
			Av.add(indexBuffer[i] * 2 + 1, dvDotResponseLocal.get(2 * i + 1));
		}
		
		response_repulse.mult(dvLocal, dvDotResponseLocal);
		dvDotResponseLocal.scale(-stepSize * stepSize);
		
		for (int i = 0; i < 6; i++) {
			Av.add(indexBuffer[i] * 2, dvDotResponseLocal.get(2 * i));
			Av.add(indexBuffer[i] * 2 + 1, dvDotResponseLocal.get(2 * i + 1));
		}
	}
	
	private static DenseVector dfdx_v = new DenseVector(12);
	
	public void applyForce(double timestep, boolean isImplicit) {
		response_viscous.mult(velocities, forces_viscous);
		response_repulse.mult(positions, forces_repulsive);
		
		if (isImplicit) {
			dfdx_v.zero();
			response_repulse.mult(velocities, dfdx_v);
			dfdx_v.scale(timestep);
			forces_repulsive.add(dfdx_v);
		}
		
		int index = 0;
		t1.A.f.x += forces_viscous.get(index++);
		t1.A.f.y += forces_viscous.get(index++);
		t1.B.f.x += forces_viscous.get(index++);
		t1.B.f.y += forces_viscous.get(index++);
		t1.C.f.x += forces_viscous.get(index++);
		t1.C.f.y += forces_viscous.get(index++);
		t2.A.f.x += forces_viscous.get(index++);
		t2.A.f.y += forces_viscous.get(index++);
		t2.B.f.x += forces_viscous.get(index++);
		t2.B.f.y += forces_viscous.get(index++);
		t2.C.f.x += forces_viscous.get(index++);
		t2.C.f.y += forces_viscous.get(index++);
		
		index = 0;
		t1.A.f.x += forces_repulsive.get(index++);
		t1.A.f.y += forces_repulsive.get(index++);
		t1.B.f.x += forces_repulsive.get(index++);
		t1.B.f.y += forces_repulsive.get(index++);
		t1.C.f.x += forces_repulsive.get(index++);
		t1.C.f.y += forces_repulsive.get(index++);
		t2.A.f.x += forces_repulsive.get(index++);
		t2.A.f.y += forces_repulsive.get(index++);
		t2.B.f.x += forces_repulsive.get(index++);
		t2.B.f.y += forces_repulsive.get(index++);
		t2.C.f.x += forces_repulsive.get(index++);
		t2.C.f.y += forces_repulsive.get(index++);
	}

}
