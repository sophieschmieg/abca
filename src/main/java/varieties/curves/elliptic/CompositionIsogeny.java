package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class CompositionIsogeny<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
	private Isogeny<T> firstIsogeny;
	private Isogeny<T> secondIsogeny;
	private ProjectiveMorphism<T> asMorphism;

	public CompositionIsogeny(Isogeny<T> firstIsogeny, Isogeny<T> secondIsogeny) {
		if (!firstIsogeny.getRange().equals(secondIsogeny.getDomain())) {
			throw new ArithmeticException("Isogenies cannot be composed!");
		}
		this.firstIsogeny = firstIsogeny;
		this.secondIsogeny = secondIsogeny;
	}
	
	@Override
	public String toString() {
		return firstIsogeny.toString() + " o " + secondIsogeny.toString();
	}
	
	@Override
	public int compareTo(Isogeny<T> o) {
		if (o instanceof CompositionIsogeny<?>) {
			CompositionIsogeny<T> other = (CompositionIsogeny<T>) o;
			int cmp = firstIsogeny.compareTo(other.firstIsogeny);
			if (cmp != 0) {
				return cmp;
			}
			return secondIsogeny.compareTo(other.secondIsogeny);
		}
		return getClass().getName().compareTo(o.getClass().getName());
	}

	@Override
	public EllipticCurve<T> getDomain() {
		return firstIsogeny.getDomain();
	}

	@Override
	public EllipticCurve<T> getRange() {
		return secondIsogeny.getRange();
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
		return secondIsogeny.evaluate(firstIsogeny.evaluate(point));
	}

	@Override
	public BigInteger getDegree() {
		return firstIsogeny.getDegree().multiply(secondIsogeny.getDegree());
	}

	@Override
	public Isogeny<T> getDual() {
		return new CompositionIsogeny<T>(secondIsogeny.getDual(), firstIsogeny.getDual());
	}

	@Override
	public ProjectiveMorphism<T> asMorphism() {
		if (asMorphism == null) {
			List<Polynomial<T>> asPolynomials = new ArrayList<>();
			PolynomialRing<T> domainRing = getDomain().asGenericProjectiveScheme().homogenousPolynomialRing();
			for (Polynomial<T> second : secondIsogeny.asMorphism().asPolynomials()) {
				asPolynomials.add(domainRing.substitute(second, firstIsogeny.asMorphism().asPolynomials()));
			}
			this.asMorphism = new ProjectiveMorphism<>(getDomain().asGenericProjectiveScheme(), getRange().asGenericProjectiveScheme(), asPolynomials);
		}
		return asMorphism;
	}
}
