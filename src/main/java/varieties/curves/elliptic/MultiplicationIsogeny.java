package varieties.curves.elliptic;

import java.math.BigInteger;

import fields.interfaces.Element;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class MultiplicationIsogeny<T extends Element<T>> implements Isogeny<T> {
	private BigInteger n;
	private EllipticCurve<T> curve;
	private ProjectiveMorphism<T> asMorphism;
	
	public MultiplicationIsogeny(BigInteger n, EllipticCurve<T> curve) {
		this.n = n;
		this.curve = curve;
	}

	@Override
	public int compareTo(Isogeny<T> o) {
		if (o instanceof MultiplicationIsogeny<?>) {
			MultiplicationIsogeny<T> other =(MultiplicationIsogeny<T>)o;
			int cmp = curve.compareTo(other.curve);
			if (cmp != 0) {
				return cmp;
			}
			return n.compareTo(other.n);
		}
		return getClass().getName().compareTo(o.getClass().getName());
	}
	
	@Override
	public String toString() {
		return "[" + n + "]";
	}

	@Override
	public EllipticCurve<T> getDomain() {
		return curve;
	}

	@Override
	public EllipticCurve<T> getRange() {
		return curve;
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
		return curve.multiply(n, point);
	}

	@Override
	public BigInteger getDegree() {
		return n.multiply(n);
	}

	@Override
	public Isogeny<T> getDual() {
		return this;
	}
	
	@Override
	public ProjectiveMorphism<T> asMorphism() {
		if (asMorphism == null) {
			asMorphism = curve.multiplicationMorphism(n.intValueExact());
		}
		return asMorphism;
	}

}
