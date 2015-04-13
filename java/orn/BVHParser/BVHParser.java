import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Stack;
import java.util.StringTokenizer;

import Jama.Matrix;

/*
 * Created on 12.11.2004
 */

public class BVHParser
{
  private static final String HIERARCHY    = "HIERARCHY";
  private static final String MOTION       = "MOTION";
  private static final String ROOT         = "ROOT";
  private static final String JOINT        = "JOINT";
  private static final String OFFSET       = "OFFSET";
  private static final String END_0        = "End";
  private static final String END_1        = "Site";
  private static final String PARENT_CLOSE = "}";
  
  private String state;
  
  private String lastToken = "";
  private Joint poppedJoint = null;
  
  private ArrayList jointsList = new ArrayList();
  private Stack jointsStack = new Stack();
  private int jointNo = 1;
  
  private ArrayList uniqueJoints = new ArrayList();
  private ArrayList leafJoints = new ArrayList();
  private int connectionMatrix[][];
  
  Joint joint = null;
  private float currentX = 0;
  private float currentY = 0;
  private float currentZ = 0;
  
  private boolean chainEnd = false;
  private boolean useEndSite = false;
  
  private ArrayList motionCoords = new ArrayList();
  
  public BVHParser( String filename, String useEndSite )
  {
    if( useEndSite.equalsIgnoreCase( "y" ) )
      this.useEndSite = true;
    
    createOrderMatrix( filename );
    
    createConnectionMatrix();
    
    printJointsToFile( filename );
  }
  
  private void createConnectionMatrix()
  {
    connectionMatrix = new int[jointNo][jointNo];
    
    Joint tempLeaf;
    for( int i = 0; i < leafJoints.size(); i++ )
      {
	tempLeaf = (Joint)leafJoints.get( i );
	createConnection( tempLeaf );
      }
  }
  
  private void createConnection( Joint iJoint )
  {
    Joint[] parents;
    Joint tempParent;
    parents = iJoint.getParents();
    
    if( parents.length > 0 )
      {
	tempParent = parents[parents.length - 1];
	System.out.println( " --> CREATING CONNECTION BETWEEN [" + iJoint.getName() + "] AND [" + tempParent.getName() + "]" );
	connectionMatrix[iJoint.getNumber() - 1][tempParent.getNumber() - 1] = 1;
	createConnection( tempParent );
      }
    
  }
  
  private void printJointsToFile( String parsedFilename )
  {
    if( parsedFilename.length() > 4 )
      {
	String lastFour = parsedFilename.substring( parsedFilename.length() - 4, parsedFilename.length() );
	
	if( lastFour.equalsIgnoreCase( ".bvh" ) )
	  {
	    parsedFilename = parsedFilename.substring( 0, parsedFilename.length() - 4 );
	  }
	
	parsedFilename += ".txt";
      }
    
    try
      {
	//PrintWriter out_joints = new PrintWriter( new BufferedWriter( new FileWriter( "JointCoordinates_" + parsedFilename ) ) );
	//PrintWriter out_jointNumbers = new PrintWriter( new BufferedWriter( new FileWriter( "JointOrder_" + parsedFilename ) ) );
	PrintWriter out_motionNumbers = new PrintWriter( new BufferedWriter( new FileWriter( "output\\Motion_" + parsedFilename ) ) );
	PrintWriter out_connectionMatrix = new PrintWriter( new BufferedWriter( new FileWriter( "output\\Connection_" + parsedFilename ) ) );
	
	//			for( int i = 0; i < jointsList.size(); i++ )
	//			{
	//				Joint joint = (Joint)jointsList.get( i );
	//				out_joints.println( joint.getX() + " " + joint.getY() + " " + joint.getZ() );
	//				out_jointNumbers.print( joint.getNumber() );
	//				if( i < jointsList.size() - 1 )
	//					out_jointNumbers.print( " " );
	//			}
	
	System.out.println( "" );
	//System.out.println( "-> WRITTEN JOINT ORDER TO [JointOrder_" + parsedFilename + "]" );
	//System.out.println( "-> WRITTEN JOINT COORDINATES TO [JointCoordinates_" + parsedFilename + "]" );
	
	if( motionCoords.size() > 0 )
	  {
	    for( int i = 0; i < motionCoords.size(); i++ )
	      {
		out_motionNumbers.println( ((StringBuffer)motionCoords.get( i )).toString() );
	      }
	    
	    System.out.println( "-> WRITTEN MOTION COORDINATES ORDER TO [Motion_" + parsedFilename + "]" );
	  }
	
	for( int i = 0; i < connectionMatrix.length; i++ )
	  {
	    for( int j = 0; j < connectionMatrix.length; j++ )
	      {
		out_connectionMatrix.print( connectionMatrix[i][j] );
		
		if( j < connectionMatrix.length - 1 )
		  out_connectionMatrix.print( " " );
	      }
	    
	    out_connectionMatrix.println( "\n" );
	  }
	
	System.out.println( "-> WRITTEN JOINT CONNECTION MATRIX TO [Connection_" + parsedFilename + "]" );
	
	//			out_jointNumbers.flush();
	//			out_jointNumbers.close();
	//			
	//			out_joints.flush();
	//			out_joints.close();
	
	out_motionNumbers.flush();
	out_motionNumbers.close();
	
	out_connectionMatrix.flush();
	out_connectionMatrix.close();
      }
    catch( IOException e )
      {
	System.out.println( "ERROR WRITING TO FILES : " + e );
      }
  }
  
  private void addToJointListAndIncrement( Joint joint )
  {
    System.out.println( "JOINT[" + joint.getName() + "] HAS [" + jointsStack.size() + "] PARENTS" );
    Joint parents[] = new Joint[jointsStack.size()];
    for( int i = 0; i < jointsStack.size(); i++ )
      {
	Joint parentJoint = (Joint)jointsStack.elementAt( i );
	parents[i] = parentJoint;
	joint.setParents( parents );
      }
    
    System.out.println( "ADDING [" + joint.getName() + "] #[" + joint.getNumber() + "] POS[" + joint.getX() + " " + joint.getY() + " " + joint.getZ() + "]" );
    jointsList.add( joint );
    uniqueJoints.add( joint );
    this.jointNo++;
  }
  
  private void addBacktrackToJointList( Joint joint )
  {
    System.out.println( "BACKTRACKING TO [" + joint.getName() + "] #[" + joint.getNumber() + "]" );
    jointsList.add( joint );		
  }
  
  private void createOrderMatrix( String filename )
  {
    // Open file
    EasyReader reader = new EasyReader( "bvh\\" + filename );
    String line = "";
    String splitted[];
    
    while( !reader.eof() )
      {
	line = ( reader.readString() ).trim();
	//System.out.println( line );
	
	if( !line.equals( "" ) )
	  {
	    splitted = splitString( line );
	    
	    if( splitted[0].equals( HIERARCHY ) )
	      {
		state = HIERARCHY;
	      }
	    else if( splitted[0].equals( MOTION ) )
	      {
		state = MOTION;	
		reader.readString();	// Skip Frames
		reader.readString();	// Skip Frame Time
	      }
	    else if( state.equals( HIERARCHY ) )
	      {
		if( splitted[0].equals( ROOT ) || splitted[0].equals( JOINT ) )
		  {
		    if( lastToken.equals( PARENT_CLOSE ) )
		      {
			if( jointsStack.size() > 0 )
			  {
			    Joint peekJoint = (Joint)jointsStack.peek();
			    addBacktrackToJointList( peekJoint );
			    
			    currentX = peekJoint.getX();
			    currentY = peekJoint.getY();
			    currentZ = peekJoint.getZ();			
			  }				
		      }
		    
		    System.out.println( "CREATING OBJECT [" + splitted[1] + "]" );
		    joint = new Joint( jointNo, splitted[1] );
		  }
		else if( splitted[0].equals( END_0 ) && splitted[1].equals( END_1 ) )
		  {
		    /* Put last joint as leaf joint */
		    leafJoints.add( joint );
		    
		    //if( useEndSite )
		    System.out.println( "CREATING OBJECT [" + END_0 + "]" );
		    joint = new Joint( jointNo, END_0 );
		    
		    chainEnd = true;
		  }
		else if( splitted[0].equals( OFFSET ) )
		  {
		    float offsetX = Float.parseFloat( splitted[1] );
		    float offsetY = Float.parseFloat( splitted[2] );	// Ath hér víxla 3 og 2
		    float offsetZ = Float.parseFloat( splitted[3] );
		    
		    currentX += offsetX;
		    currentY += offsetY;
		    currentZ += offsetZ;
		    
		    // Read offset values
		    joint.setOffsetX( offsetX );
		    joint.setOffsetY( offsetY );
		    joint.setOffsetZ( offsetZ );
		    joint.setX( currentX );
		    joint.setY( currentY );
		    joint.setZ( currentZ );
		    
		    // Submit Joint to collection
		    if( !joint.getName().equals( END_0 ) || ( useEndSite && joint.getName().equals( END_0 ) ) )
		      addToJointListAndIncrement( joint );
		    
		    jointsStack.push( joint );
		  }
		else if( splitted[0].equals( PARENT_CLOSE ) )
		  {
		    // Pop current parent if is end of chain otherwise it will backtrack to itself
		    if( chainEnd )
		      {
			jointsStack.pop();
			chainEnd = false;
		      }
		    else
		      {
			if( jointsStack.size() > 0 )
			  {
			    poppedJoint = (Joint)jointsStack.pop();
			    addBacktrackToJointList( poppedJoint );
			    
			    if( jointsStack.size() > 0 )
			      {
				currentX = ( (Joint)jointsStack.peek() ).getX();
				currentY = ( (Joint)jointsStack.peek() ).getY();
				currentZ = ( (Joint)jointsStack.peek() ).getZ();
			      }
			  }
		      }				
		  }
		
		lastToken = splitted[0];
	      }
	    else if( state.equals( MOTION ) )
	      {
		ArrayList frameVertices = new ArrayList();
		int currentJointNo = 0;
		
		// Ignore 0-2 to make character stay in place
		int pos = 3;
		
		// Currently only supports order ZXY
		for( int i = 3; i < splitted.length; )
		  {
		    double rotZ = Float.parseFloat( splitted[i] );
		    rotZ = degreesToRad( rotZ );
		    i++;
		    double rotX = Float.parseFloat( splitted[i] );	// Swap X and Y
		    rotX = degreesToRad( rotX );
		    i++;
		    double rotY = Float.parseFloat( splitted[i] );
		    rotY = degreesToRad( rotY );
		    i++;
		    
		    double[][] coordsX = { { 1, 0, 0, 0 }, { 0, Math.cos(rotX), -Math.sin(rotX), 0 }, { 0, Math.sin(rotX), Math.cos(rotX), 0 }, { 0, 0, 0, 1 } };
		    double[][] coordsY = { { Math.cos(rotY), 0, Math.sin(rotY), 0 }, { 0, 1, 0, 0 }, { -Math.sin(rotY), 0, Math.cos(rotY), 0 }, { 0, 0, 0, 1 } };
		    double[][] coordsZ = { { Math.cos(rotZ), -Math.sin(rotZ), 0, 0 }, { Math.sin(rotZ), Math.cos(rotZ), 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
		    double[][] translate = { { 0, 0, 0, (double)((Joint)uniqueJoints.get( currentJointNo )).getOffsetX() },
					     { 0, 0, 0, (double)((Joint)uniqueJoints.get( currentJointNo )).getOffsetY() },
					     { 0, 0, 0, (double)((Joint)uniqueJoints.get( currentJointNo )).getOffsetZ() },
					     { 0, 0, 0, 0 } };
		    
		    Matrix matrixX = new Matrix( coordsX );
		    Matrix matrixY = new Matrix( coordsY );
		    Matrix matrixZ = new Matrix( coordsZ );
		    Matrix matrixTranslate = new Matrix( translate );
		    
		    Matrix matrixM = matrixZ;
		    matrixM = matrixM.times( matrixX );
		    matrixM = matrixM.times( matrixY );
		    matrixM = matrixM.plus( matrixTranslate );
		    
		    frameVertices.add( matrixM );
		    
		    currentJointNo++;
		  }
		
		for( int i = 0; i < frameVertices.size(); i++ )
		  {
		    Matrix m;
		    double mult[] = { 0, 0, 0, 1 };
		    Matrix rightMost = new Matrix( mult, 1 );
		    rightMost = rightMost.transpose();
		    
		    Joint tempJoint = (Joint)uniqueJoints.get( i );
		    Joint parents[] = tempJoint.getParents();
		    
		    if( parents.length > 0 )
		      {
			
			m = (Matrix)frameVertices.get( parents[0].getNumber() - 1 );
			
			for( int j = 1; j < parents.length; j++ )
			  {
			    m = m.times( (Matrix)frameVertices.get( parents[j].getNumber() - 1 ) );
			  }
			
			m = m.times( (Matrix)frameVertices.get( i ) );
		      }
		    else
		      {
			m = (Matrix)frameVertices.get( i );
		      }
		    
		    m = m.times( rightMost );
		    
		    double results[][] = m.getArray();
		    tempJoint.setMotionX( (float)results[0][0] );
		    tempJoint.setMotionY( (float)results[1][0] );
		    tempJoint.setMotionZ( (float)results[2][0] );
		  }
		
		addToList();
	      }
	    else
	      {
		System.out.println( "ERROR IN FILE!" );
	      }
	  }
      }
  }
  
  private void addToList()
  {
    StringBuffer cJoints = new StringBuffer();
    
    Joint tempJoint;
    for( int i = 0; i < uniqueJoints.size(); i++ )
      {
	tempJoint = (Joint)uniqueJoints.get( i );
	cJoints.append( tempJoint.getMotionX() + " " );
      }
    
    for( int i = 0; i < uniqueJoints.size(); i++ )
      {
	tempJoint = (Joint)uniqueJoints.get( i );
	cJoints.append( tempJoint.getMotionY() + " " );
      }
    
    for( int i = 0; i < uniqueJoints.size(); i++ )
      {
	tempJoint = (Joint)uniqueJoints.get( i );
	cJoints.append( tempJoint.getMotionZ() );
	if( i < uniqueJoints.size() - 1 )
	  cJoints.append( " " );
      }
    
    motionCoords.add( cJoints );
  }
  
  private double degreesToRad( double degrees )
  {
    return ( degrees * Math.PI ) / 180.0;
  }
  
  private String[] splitString( String line )
  {
    StringTokenizer tokenizer = new StringTokenizer( line );
    String splitted[] = new String[tokenizer.countTokens()];
    
    for( int i = 0; i < splitted.length; i++ )
      {
	splitted[i] = tokenizer.nextToken();
      }
    
    return splitted;
  }
  
  public static void main( String args[] )
  {
    String useEnd;
    if( args.length > 1 )
      useEnd = args[1];
    else
      useEnd = "n";
    
    
    BVHParser parser = new BVHParser( args[0], useEnd );
  }
}
