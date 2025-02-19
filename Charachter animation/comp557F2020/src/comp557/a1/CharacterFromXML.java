package comp557.a1;
 		  	  				   
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.Scanner;

import javax.vecmath.Tuple2d;
import javax.vecmath.Tuple3d;
import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

//**** Marwen Kraiem 260955202 ****//


/**
 * Loads an articulated character hierarchy from an XML file. 
 */
public class CharacterFromXML {

	public static GraphNode load( String filename ) {
		try {
			InputStream inputStream = new FileInputStream(new File(filename));
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder builder = factory.newDocumentBuilder();
			Document document = builder.parse(inputStream);
			return createScene( null, document.getDocumentElement() ); // we don't check the name of the document elemnet
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException("Failed to load simulation input file.", e);
		}
	}
	
	/**
	 * Load a subtree from a XML node.
	 * Returns the root on the call where the parent is null, but otherwise
	 * all children are added as they are created and all other deeper recursive
	 * calls will return null.
	 */
	public static GraphNode createScene( GraphNode parent, Node dataNode ) {
        NodeList nodeList = dataNode.getChildNodes();
        for ( int i = 0; i < nodeList.getLength(); i++ ) {
            Node n = nodeList.item(i);
            // skip all text, just process the ELEMENT_NODEs
            if ( n.getNodeType() != Node.ELEMENT_NODE ) continue;
            String nodeName = n.getNodeName();
            GraphNode node = null;
            if ( nodeName.equalsIgnoreCase( "node" ) ) {
            	node = CharacterFromXML.createJoint( n );
            } else if ( nodeName.equalsIgnoreCase( "geom" ) ) {        		
        		node = CharacterFromXML.createGeom( n ) ;            
            } else {
            	System.err.println("Unknown node " + nodeName );
            	continue;
            }
            if ( node == null ) continue;
            // recurse to load any children of this node
            createScene( node, n );
            if ( parent == null ) {
            	// if no parent, we can only have one root... ignore other nodes at root level
            	return node;
            } else {
            	parent.add( node );
            }
        }
        return null;
	}
	
	/**​‌​​​‌‌​​​‌‌​​​‌​​‌‌‌​​‌
	 * Create a joint
	 * 
	 * TODO: Objective 8: XML, Adapt commented code in createJoint() to create your joint nodes when loading from xml
	 */
	public static GraphNode createJoint( Node dataNode ) {
		String type = dataNode.getAttributes().getNamedItem("type").getNodeValue();
		String name = dataNode.getAttributes().getNamedItem("name").getNodeValue();
		Tuple3d t;
		Tuple2d t1;
		Tuple2d t2;
		Tuple2d t3;

		if ( type.equals("freejoint") ) {
			FreeJoint joint = new FreeJoint( name );
			return joint;
		} else if ( type.equals("spherical") ) {
			// position is optional (ignored if missing) but should probably be a required attribute!​‌​​​‌‌​​​‌‌​​​‌​​‌‌‌​​‌
			// Could add optional attributes for limits (to all joints)	
			
			//angle limits
			t1 = getTuple2dAttr(dataNode,"MinMaxRx");
			if(t1 == null) t1.set(0,0);
			t2 = getTuple2dAttr(dataNode,"MinMaxRy");
			if(t2 == null) t2.set(0,0);
			t3 = getTuple2dAttr(dataNode,"MinMaxRz");	
			if(t3 == null) t3.set(0,0);

			SphericalJoint joint = new SphericalJoint( name, 0, 0, 0, t1.x, t1.y, t2.x, t2.y, t3.x, t3.y);
			if ( (t=getTuple3dAttr(dataNode,"position")) != null ) joint.setPosition( t );

			return joint;
			
		} else if ( type.equals("rotary") ) {
			// position and axis are required... passing null to set methods
			// likely to cause an execption (perhaps OK)
			t1 = getTuple2dAttr(dataNode,"MinMaxR");
			if(t1 == null) t1.set(0,0);
			RotaryJoint joint = new RotaryJoint( name, 0, 0, 0, 0, 0, 0, t1.x, t1.y );
			joint.setPosition( getTuple3dAttr(dataNode,"position") );
			joint.setAxis( getTuple3dAttr(dataNode,"axis") );
			return joint;
			
		} else {
			System.err.println("Unknown type " + type );
		}
		return null;
	}

	/**
	 * Creates a geometry DAG node 
	 * 
	 * TODO: Objective 5: Adapt commented code in createGeom to create your geometry nodes when loading from xml
	 */
	public static GraphNode createGeom( Node dataNode ) {
		String type = dataNode.getAttributes().getNamedItem("type").getNodeValue();
		String name = dataNode.getAttributes().getNamedItem("name").getNodeValue();
		Tuple3d t;
		if ( type.equals("box" ) ) {
			BodyBox geom = new BodyBox( name, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0f, 0.0f, 0.0 );
			if ( (t=getTuple3dAttr(dataNode,"center")) != null ) geom.setCentre( t );
			if ( (t=getTuple3dAttr(dataNode,"scale")) != null ) geom.setScale( t );
			if ( (t=getTuple3dAttr(dataNode,"color")) != null ) geom.setColor( t );
			return geom;
		} else if ( type.equals( "sphere" )) {
			BodySphere geom = new BodySphere( name, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0 );				
			if ( (t=getTuple3dAttr(dataNode,"center")) != null ) geom.setCentre( t );
			if ( (t=getTuple3dAttr(dataNode,"scale")) != null ) geom.setScale( t );
			if ( (t=getTuple3dAttr(dataNode,"color")) != null ) geom.setColor( t );
			return geom;	
		} else {
			System.err.println("unknown type " + type );
		}
		return null;		
	}
	
	/**
	 * Loads tuple3d attributes of the given name from the given node.
	 * @param dataNode
	 * @param attrName
	 * @return null if attribute not present
	 */
	public static Tuple3d getTuple3dAttr( Node dataNode, String attrName ) {
		Node attr = dataNode.getAttributes().getNamedItem( attrName);
		Vector3d tuple = null;
		if ( attr != null ) {
			Scanner s = new Scanner( attr.getNodeValue() );
			tuple = new Vector3d( s.nextDouble(), s.nextDouble(), s.nextDouble() );			
			s.close();
		}
		return tuple;
	}
	
	/**
	 * Loads tuple2d attributes of the given name from the given node.
	 * @param dataNode
	 * @param attrName
	 * @return null if attribute not present
	 */
	public static Tuple2d getTuple2dAttr( Node dataNode, String attrName ) {
		Node attr = dataNode.getAttributes().getNamedItem( attrName);
		Vector2d tuple = null;
		if ( attr != null ) {
			Scanner s = new Scanner( attr.getNodeValue() );
			tuple = new Vector2d( s.nextDouble(), s.nextDouble());			
			s.close();
		}
		return tuple;
	}
	
	

}