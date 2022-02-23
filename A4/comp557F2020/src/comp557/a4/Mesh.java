/*
 * Marwen Kraiem 260955202
 */

package comp557.a4;

import java.util.HashMap;
import java.util.Map;

import javax.vecmath.Vector3d;

public class Mesh extends Intersectable {
	
	/** Static map storing all meshes by name */
	public static Map<String,Mesh> meshMap = new HashMap<String,Mesh>();
	
	/**  Name for this mesh, to allow re-use of a polygon soup across Mesh objects */
	public String name = "";
	
	/**
	 * The polygon soup.
	 */
	public PolygonSoup soup;

	public Mesh() {
		super();
		this.soup = null;
	}			
		
	@Override
	public void intersect(Ray ray, IntersectResult result) {
		
		// TODO: Objective 9: ray triangle intersection for meshes
		Vector3d ab = new Vector3d();
		Vector3d bc = new Vector3d();
		Vector3d ca = new Vector3d();

		Vector3d ax = new Vector3d();
		Vector3d n = new Vector3d();
		
		Vector3d p = new Vector3d(ray.eyePoint);
		
		double t;
		
		for(int[] face:soup.faceList) {
			
			ab.sub(soup.vertexList.get(face[1]).p,soup.vertexList.get(face[0]).p);
			bc.sub(soup.vertexList.get(face[2]).p,soup.vertexList.get(face[1]).p);
			ca.sub(soup.vertexList.get(face[0]).p,soup.vertexList.get(face[2]).p);
			
			n.cross(ab, bc);
			p.sub(soup.vertexList.get(face[0]).p, ray.eyePoint);
			if(n.length()<1e-9)  n.cross(bc, ca);
			if(n.length()<1e-9) n.cross(ca, ab);
			if(n.length()>1e-9)	n.normalize(); else continue;
			
			
			if(Math.abs(ray.viewDirection.dot(n))>1e-9) {
			    t = p.dot(n)/ray.viewDirection.dot(n); 	
			    
			    ax.scaleAdd(t,ray.viewDirection,ray.eyePoint);
			    ax.sub(soup.vertexList.get(face[0]).p);	
			    ax.cross(ab, ax);
			    
			    if (ax.dot(n)>0){
				    
			    	ax.scaleAdd(t,ray.viewDirection,ray.eyePoint);
				    ax.sub(soup.vertexList.get(face[1]).p);	
				    ax.cross(bc, ax);
				    
			    	if (ax.dot(n)>0) {
					    
			    		ax.scaleAdd(t,ray.viewDirection,ray.eyePoint);
					    ax.sub(soup.vertexList.get(face[2]).p);	
					    ax.cross(ca, ax);
			    		
			    		if(ax.dot(n)>0) {
			    			
			    			if (t > 1e-9 && t < result.t) { 
									result.t = t;
									result.p.scaleAdd(t,ray.viewDirection,ray.eyePoint);	
									result.n.set(n); 
									result.material = material;
							}
			    		
			    		}
			    	
			    	}
			    	
			    }
			}
		}
		
	}

}
