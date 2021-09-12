package varieties;

import fields.helper.CoordinateRing;
import fields.interfaces.Element;
import varieties.affine.AffineCover;
import varieties.affine.AffineScheme;

public abstract class AbstractScheme<T extends Element<T>, S extends Element<S>> implements Scheme<T, S> {

	@Override
	public boolean isConnected() {
		AffineCover<T> affineCover = getAffineCover();
		for (int i = 0; i < affineCover.getCover().size(); i++){
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

}
