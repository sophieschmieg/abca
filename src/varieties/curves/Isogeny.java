package varieties.curves;

import java.math.BigInteger;

import fields.Element;
import varieties.ProjectivePoint;

public interface Isogeny<T extends Element> {
	public EllipticCurve<T> getDomain();
	public EllipticCurve<T> getRange();
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point);
	public BigInteger getDegree();
	public Isogeny<T> getDual();
}
