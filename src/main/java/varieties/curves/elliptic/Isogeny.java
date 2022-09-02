package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.MathMap;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public interface Isogeny<T extends Element<T>>
		extends MathMap<ProjectivePoint<T>, ProjectivePoint<T>>, Element<Isogeny<T>> {
	public BigInteger getDegree();

	public EllipticCurve<T> getDomain();

	public EllipticCurve<T> getRange();

	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point);

	public ProjectiveMorphism<T> asMorphism();

	public Isogeny<T> getDual();
	
	public List<ProjectivePoint<T>> kernelGenerators();
}
