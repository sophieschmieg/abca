package fields.quaternions;

import fields.interfaces.Algebra;
import fields.interfaces.DedekindRing;
import fields.interfaces.Element;
import fields.quaternions.AbstractQuaternions.Quaternion;

public interface QuaternionOrder<T extends Element<T>, I extends Element<I>, R extends Element<R>> extends Algebra<I, Quaternion<T>> {
	boolean isElement(Quaternion<T> t);
	
	I reducedDiscriminant();
	boolean isMaximal();
	
	DedekindRing<I, T, R> getRing();
	Quaternions<T> getQuaternions();
}
