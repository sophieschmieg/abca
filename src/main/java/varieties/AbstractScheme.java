package varieties;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;

import fields.interfaces.Element;
import fields.polynomials.CoordinateRing;
import varieties.affine.AffineCover;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;

public abstract class AbstractScheme<T extends Element<T>, S extends Element<S>> implements Scheme<T, S> {

	@Override
	public boolean isConnected() {
		AffineCover<T> affineCover = getAffineCover();
		for (int i = 0; i < affineCover.getCover().size(); i++) {
			for (int j = 0; j < affineCover.getCover().size(); j++) {
				CoordinateRing<T> coordinateRing = affineCover.getIntersection(i, j).getCoordinateRing();
				if (coordinateRing.getIdeal().contains(coordinateRing.getPolynomialRing().one())) {
					return false;
				}
			}
		}
		return true;
	}

	@Override
	public boolean isIntegral() {
		if (!isConnected()) {
			return false;
		}
		for (AffineScheme<T> affine : getAffineCover().getCover()) {
			if (!affine.isIntegral()) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean isReduced() {
		for (AffineScheme<T> affine : getAffineCover().getCover()) {
			if (!affine.isReduced()) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean isIrreducible() {
		if (!isConnected()) {
			return false;
		}
		for (AffineScheme<T> affine : getAffineCover().getCover()) {
			if (!affine.isIrreducible()) {
				return false;
			}
		}
		return true;
	}

	@Override
	public final int dimension() {
		return getAffineCover().getCover().get(0).getCoordinateRing().krullDimension();
	}

	@Override
	public int degree() {
		return getAffineCover().getCover().get(0).getCoordinateRing().degree();
	}

	@Override
	public List<S> singularPoints() {
		Set<S> candidates = new TreeSet<>();
		for (int i = 0; i < getAffineCover().getCover().size(); i++) {
			AffineScheme<T> affine = getAffineCover().getCover().get(i);
			for (AffinePoint<T> point : affine.singularPoints()) {
				candidates.add(fromAffinePoint(point, i));
			}
		}
		List<S> result = new ArrayList<>();
		result.addAll(candidates);
		return result;
	}

	@Override
	public Optional<? extends Morphism<T, S, S>> singularLocus() {
return null;
		//		List<Morphism<T, AffinePoint<T>, S>> result = new ArrayList<>();
//		for (int i = 0; i < getAffineCover().getCover().size(); i++) {
//			AffineScheme<T> affine = getAffineCover().getCover().get(i);
//			for (Morphism<T, AffinePoint<T>, AffinePoint<T>> singular : affine.singularLocus()) {
//				result.add(new CompositionMorphism<>(singular, embedding(i)));
//			}
//		}
//		return result;
	}
}
