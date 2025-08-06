package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.List;

import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field.Extension;
import fields.interfaces.FieldExtension;
import fields.interfaces.MathMap;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public interface Isogeny<T extends Element<T>>
		extends MathMap<ProjectivePoint<T>, ProjectivePoint<T>>, Element<Isogeny<T>> {
	public BigInteger getDegree();

	public EllipticCurve<T> getDomain();

	public EllipticCurve<T> getRange();

	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point);
	
	public <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> ProjectivePoint<E> evaluateOverExtension(ProjectivePoint<E>point, Extension<T, B, E, FE> extension);

	public ProjectiveMorphism<T> asMorphism();

	public Isogeny<T> getDual();
	
	public List<ProjectivePoint<T>> kernelGenerators();
}
