/*
 * Marwen Kraiem 260955202
 */

package comp557.a4;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.vecmath.Color3f;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * Simple scene loader based on XML file format.
 */
public class Scene {
    
    /** List of surfaces in the scene */
    public List<Intersectable> surfaceList = new ArrayList<Intersectable>();
	
	/** All scene lights */
	public Map<String,Light> lights = new HashMap<String,Light>();

    /** Contains information about how to render the scene */
    public Render render;
    
    /** The ambient light colour */
    public Color3f ambient = new Color3f();

    
    /** Maximum recursion depth **/
    public  int MaxDepth = 1;
    
    /** 
     * Default constructor.
     */
    public Scene() {
    	this.render = new Render();
    }
   
    /**
     * renders the scene
     */
    public void render(boolean showPanel) {
 
        Camera cam = render.camera; 
        int w = cam.imageSize.width;
        int h = cam.imageSize.height;
        
               
        render.init(w, h, showPanel);
        
        Color3f c = new Color3f();
        Color3f tmpc = new Color3f();
        
        Ray ray = new Ray();
        IntersectResult result = new IntersectResult();
        
        double[] offset = new double[2];
        int n = (int) Math.sqrt(render.samples); 
        int m = (int)render.samples/n;
        
        // rendering loops 
        for ( int j = 0; j < h && !render.isDone(); j++ ) {
            for ( int i = 0; i < w && !render.isDone(); i++ ) {
            	
            	// initialize color	            	
            	c.set(0,0,0);
            	for (int e = 0; e<cam.eyeSamples; e++) {
            		tmpc.set(0,0,0);
	            	for (int p = 0; p < n; p++) {
	            		offset[0] = render.jitter? (p+Math.random()-0.5)/(float)n: (p+0.5)/(float)n; 
	            		for (int q = 0; q <m ; q++) {
	            			offset[1] = render.jitter? (q+Math.random()-0.5)/(float)m: (q+0.5)/(float)m;
	            			
		        		   // TODO: Objective 1: generate a ray (use the generateRay method) 
	            			ray = new Ray();
	            			generateRay( i, j,offset, cam,ray);
				        		   
			               // TODO: Objective 2: test for intersection with scene surfaces
	            		   result = new IntersectResult();
				           for (Intersectable surface : surfaceList)  surface.intersect(ray, result);			           
			               				           
				           // TODO: Objective 3: compute the shaded result for the intersection point (perhaps requiring shadow rays)
				           if (result.t < Double.POSITIVE_INFINITY ) {
				        	    tmpc.add(generateColor(result, ray,1,MaxDepth)); 
				            }      
				           else 
				        	   tmpc.add(render.bgcolor);			           
	                   }
		           }
	            	tmpc.scale(1.0f /(float)(n*m));
	            	c.add(tmpc);            	

            	}
            	c.scale(1.0f/(float)(cam.eyeSamples));
            
            	// clamp color between 0 and 1
	           c.clamp(0,1);
	            
	           int r = (int)(255*c.x);
	           int g = (int)(255*c.y);
               int b = (int)(255*c.z);
               int a = 255;
               int argb = (a<<24 | r<<16 | g<<8 | b);    
           
                // update the render image	
                render.setPixel(i, j, argb);     

            }
       }
        
        // save the final render image
        render.save();
        
        // wait for render viewer to close
        render.waitDone();
        
    }
    
    
    /***
     *  compute color from a ray and intersection results
     * @param result
     * @param ray
     * @param n : refraction index of the space
     * @param depth  : recursion depth
     */
    private Color3f generateColor(IntersectResult result , Ray ray, double n,int depth) {
	         Color3f tmp = new Color3f();
	         Color3f color = new Color3f();
	         Color3f c = new Color3f();

	         double diffuse;
	         double specular;	
	         // view direction
	         Vector3d v = new Vector3d();
	         
	 		// half-direction
	 		Vector3d half = new Vector3d();
	    		    				
	 		Ray shadowRay = new Ray();
	    
	 		//ambient color		 	  	
	    	color.set(result.material.diffuse.x*ambient.x,result.material.diffuse.y*ambient.y,result.material.diffuse.z*ambient.z);
	    	
	    	for( Light light :lights.values()) {
	    		c.set(0,0,0);
			    for (float i = 0; i<light.samples;i++){ 
						// generate a shadow ray   
	    				shadowRay = new Ray();
	    				shadowRay.viewDirection.set(light.from);
	    				// light area
	    				if(light.type.contains("area")) {
	    					shadowRay.viewDirection.scaleAdd(Math.random(), light.u,shadowRay.viewDirection);
	    					shadowRay.viewDirection.scaleAdd(Math.random(), light.v, shadowRay.viewDirection);
	    				}
	    				shadowRay.viewDirection.sub(result.p);
	    				shadowRay.eyePoint.set(result.p);
	    				
			    		if(!inShadow( result, light, surfaceList ,new IntersectResult(), shadowRay)) {			    						    	

			    			// light direction				    			
			    			shadowRay.viewDirection.normalize();
			    						
			    			// compute the view direction    			
			    			v.scale(-1, ray.viewDirection);
			    			v.normalize();
			    			
			        		// halfVector direction
			        		half.add(v,shadowRay.viewDirection);
			        		half.normalize();
			        		
			        		diffuse = Math.max(0.0f, result.n.dot(shadowRay.viewDirection));
			        		specular = Math.max(0.0f, result.n.dot(half));
			        			
			        		if (diffuse == 0.0)
			        			specular = 0.0;
			        		else {
			            		specular = Math.pow(specular, result.material.shinyness); // sharpen the highlight
			            		
			            		// scattered Light 	
			            		tmp.set(result.material.diffuse.x*light.color.x, result.material.diffuse.y*light.color.y, result.material.diffuse.z*light.color.z);
			            		tmp.scale((float) diffuse);
			            		tmp.scale((float) light.power); 

			            		c.add(tmp);
			
			            		// reflected Light
			            		tmp.set(result.material.specular.x*light.color.x, result.material.specular.y*light.color.y, result.material.specular.z*light.color.z); ;
			            		tmp.scale((float) specular);	
			            		tmp.scale((float) light.power); 

			            		
			            		c.add(tmp);
	
			        			}


			    		}

			    }
			    c.scale(1.0f/(float)light.samples);
			    color.add(c);
	    	}
			    		
	    	if(depth>0) {
			if (result.material.mirror.w != 0) {
				Ray reflectedRay = new Ray();
			    IntersectResult reflectedResult = new IntersectResult();
			            			
			    reflectedRay.eyePoint.set(result.p);
			    reflectedRay.viewDirection.scaleAdd(-2*result.n.dot(ray.viewDirection), result.n, ray.viewDirection);		            			
			    for(Intersectable surface: surfaceList) surface.intersect(reflectedRay, reflectedResult); 
			    if (reflectedResult.t < Double.POSITIVE_INFINITY ) {
		        	tmp.set(generateColor(reflectedResult,reflectedRay,result.material.n,depth)); 
			        tmp.set(result.material.mirror.x*tmp.x,result.material.mirror.y*tmp.y,result.material.mirror.z*tmp.z);
			        color.add(tmp);
			            				
			    } else color.add(render.bgcolor);  
			}
				                       
			 if (result.material.refractive.w != 0 ) {
			 	Ray refractedRay = new Ray();
		        IntersectResult refractedResult = new IntersectResult();
		        refractedRay.eyePoint.set(result.p);
	
		        double cosi = ray.viewDirection.dot(result.n);
		        Vector3d Normal = new Vector3d();
		        Normal.set(result.n);
		        cosi =  Math.max(-1, Math.min(1, cosi));
		        double n1;
		        if (cosi<0) {
		        	cosi=-cosi; 
		        	n1 = n/result.material.n;
		        	}else {
		        		Normal.scale(-1);
		        		n1 = result.material.n/n;
		        	};
		        double s = 1-n1*n1*(1-cosi*cosi);
			   
		        if (s>=0) {
		        refractedRay.viewDirection.scale(-Math.sqrt(s), Normal);
			    refractedRay.viewDirection.scaleAdd(n1, ray.viewDirection, refractedRay.viewDirection);		            			
			    refractedRay.viewDirection.scaleAdd(n1*cosi, Normal, refractedRay.viewDirection);		            			
			    
			    
			    for(Intersectable surface: surfaceList) surface.intersect(refractedRay, refractedResult); 
			    if (refractedResult.t < Double.POSITIVE_INFINITY) {
			        	tmp.set(generateColor(refractedResult,refractedRay,result.material.n,depth)); ;
			            tmp.set(result.material.refractive.x*tmp.x,result.material.refractive.y*tmp.y,result.material.refractive.z*tmp.z);
			            color.add(tmp);
			            				
			        } else color.add(render.bgcolor);
			            						
		        }
		            					            
			 }
			 depth--;
	    	}		
			    
			return color;    	
    }
        
    /**
     * Generate a ray through pixel (i,j).
     * 
     * @param i The pixel row.
     * @param j The pixel column.
     * @param offset The offset from the center of the pixel, in the range [-0.5,+0.5] for each coordinate. 
     * @param cam The camera.
     * @param ray Contains the generated ray.
     */
	public static void generateRay(final int i, final int j, final double[] offset, final Camera cam, Ray ray) {
		
		// TODO: Objective 1: generate rays given the provided parameters
		
		double ratio = (double)cam.imageSize.width/(double)cam.imageSize.height;	
		double d  = cam.from.distance(cam.to);
		double t = Math.tan((Math.PI/180*cam.fovy)/2 )* d;;
		double b = -t ;
		double r = t*ratio;
    	double l = -r;
		
		double _u = l + (r-l)*(i + offset[0])/cam.imageSize.width;		
		double _v = -b - (t-b)*(j + offset[1])/cam.imageSize.height;	
		
		Vector3d u = new Vector3d();
	    Vector3d v = new Vector3d();
		Vector3d w = new Vector3d();
	    
		w.sub( cam.from,cam.to);
	    w.normalize();
    	u.cross(cam.up,w);
    	u.normalize();
    	v.cross(w, u);
    	    	
    	// viewDirection  
    	Vector3d viewDirection = new Vector3d();
    	
    	viewDirection.scale(_u, u);
		viewDirection.scaleAdd(_v, v, viewDirection);
		viewDirection.scaleAdd(-d, w, viewDirection);
		/**
		 * if we have camera focus (positive), we change the ray direction to be defined from the new shifted origin to the focal point. The focal point is the ray direction multiplied by the focus+camera.origin. 
		 * ref: https://medium.com/@elope139/depth-of-field-in-path-tracing-e61180417027#:~:text=Implementing%20depth%20of%20field%20in,out%20of%20focus%20will%20appear. and 13.4.3 FCG		 
		 **/
		
		if (cam.focus>0){
		  	Point3d eyeOffset = new Point3d(); 
		  	eyeOffset.set(Math.random()-0.5,Math.random()-0.5,Math.random()-0.5);// random vector with values between -0.5 and 0.5    		  	
		  	eyeOffset.scale(cam.Blur);
		  	eyeOffset.add(cam.from);	
		  	
			viewDirection.normalize();
			viewDirection.scale(cam.focus);
			viewDirection.sub(eyeOffset);
			viewDirection.add(cam.from);

			ray.set(eyeOffset, viewDirection);	

		}
		
		else {
			ray.set(cam.from,viewDirection );
			}
		
	}


	/**
	 * Shoot a shadow ray in the scene and get the result.
	 * 
	 * @param result Intersection result from ray tracing. 
	 * @param light The light to check for visibility.
	 * @param root The scene node.
	 * @param shadowResult Contains the result of a shadow ray test.
	 * @param shadowRay Contains the shadow ray used to test for visibility.
	 * 
	 * @return True if a point is in shadow, false otherwise. 
	 */
	public static boolean inShadow(final IntersectResult result, final Light light, final List<Intersectable> surfaces, IntersectResult shadowResult, Ray shadowRay) {
		
		// TODO: Objective 5: check for shadows and use it in your lighting computation
		
		
		// shadow ray intersections
		for(Intersectable surface: surfaces)  surface.intersect(shadowRay, shadowResult);
		
 		if(shadowResult.t > 0 && shadowResult.t < 1 ) {
			return true;
		}
		 return false;
	} 

}
