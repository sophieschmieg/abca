package varieties.curves;

import java.math.BigInteger;

import fields.interfaces.Element;
import fields.interfaces.MathMap;
import varieties.ProjectiveMorphism;
import varieties.ProjectivePoint;

public interface Isogeny<T extends Element<T>> extends MathMap<ProjectivePoint<T>,ProjectivePoint<T>> {
	public EllipticCurve<T> getDomain();
	public EllipticCurve<T> getRange();
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point);
	public ProjectiveMorphism<T> asMorphism();
	public BigInteger getDegree();
	public Isogeny<T> getDual();
}
