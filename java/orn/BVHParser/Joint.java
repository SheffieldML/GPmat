/*
 * Created on 12.11.2004
 */

public class Joint
{
	private int number;
	private String name;
	
	private float x;
	private float y;
	private float z;
	
	private float offsetX;
	private float offsetY;
	private float offsetZ;
	
	private float motionX;
	private float motionY;
	private float motionZ;
	
	private Joint parents[] = {};
	
	public Joint( int number )
	{
		setNumber( number );
	}
	
	public Joint( int number, String name )
	{
		setNumber( number );
		setName( name );
	}

	/**
	 * @return
	 */
	public int getNumber()
	{
		return number;
	}

	/**
	 * @return
	 */
	public float getOffsetX()
	{
		return offsetX;
	}

	/**
	 * @return
	 */
	public float getOffsetY()
	{
		return offsetY;
	}

	/**
	 * @return
	 */
	public float getOffsetZ()
	{
		return offsetZ;
	}

	/**
	 * @return
	 */
	public float getX()
	{
		return x;
	}

	/**
	 * @return
	 */
	public float getY()
	{
		return y;
	}

	/**
	 * @return
	 */
	public float getZ()
	{
		return z;
	}

	/**
	 * @param i
	 */
	public void setNumber(int i)
	{
		number = i;
	}

	/**
	 * @param f
	 */
	public void setOffsetX(float f)
	{
		offsetX = f;
	}

	/**
	 * @param f
	 */
	public void setOffsetY(float f)
	{
		offsetY = f;
	}

	/**
	 * @param f
	 */
	public void setOffsetZ(float f)
	{
		offsetZ = f;
	}

	/**
	 * @param f
	 */
	public void setX(float f)
	{
		x = f;
	}

	/**
	 * @param f
	 */
	public void setY(float f)
	{
		y = f;
	}

	/**
	 * @param f
	 */
	public void setZ(float f)
	{
		z = f;
	}

	/**
	 * @return
	 */
	public String getName()
	{
		return name;
	}

	/**
	 * @param string
	 */
	public void setName(String string)
	{
		name = string;
	}

	/**
	 * @return
	 */
	public Joint[] getParents()
	{
		return parents;
	}

	/**
	 * @param strings
	 */
	public void setParents(Joint[] i)
	{
		parents = i;
	}

	/**
	 * @return
	 */
	public float getMotionX()
	{
		return motionX;
	}

	/**
	 * @return
	 */
	public float getMotionY()
	{
		return motionY;
	}

	/**
	 * @return
	 */
	public float getMotionZ()
	{
		return motionZ;
	}

	/**
	 * @param f
	 */
	public void setMotionX(float f)
	{
		motionX = f;
	}

	/**
	 * @param f
	 */
	public void setMotionY(float f)
	{
		motionY = f;
	}

	/**
	 * @param f
	 */
	public void setMotionZ(float f)
	{
		motionZ = f;
	}

}
