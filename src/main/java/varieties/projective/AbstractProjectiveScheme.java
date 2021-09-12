package varieties.projective;

import java.util.List;

import fields.interfaces.Element;
import varieties.AbstractScheme;
import varieties.Morphism;
import varieties.SpectrumOfField.SingletonPoint;
import varieties.affine.AffineCover;
import varieties.affine.AffinePoint;

public abstract class AbstractProjectiveScheme<T extends Element<T>> extends AbstractScheme<T, ProjectivePoint<T>> implements ProjectiveVarietyInterface<T> {
	private ProjectiveScheme<T> asProjectiveVariety;

	public AbstractProjectiveScheme(ProjectiveScheme<T> asProjectiveVariety) {
		this.asProjectiveVariety = asProjectiveVariety;
	}
	
	@Override
	public Exactness exactness() {
		return asProjectiveVariety.exactness();
	}

	@Override
	public ProjectiveScheme<T> asProjectiveVariety() {
		return asProjectiveVariety;
	}

	@Override
	public AffineCover<T> getAffineCover() {
		return asProjectiveVariety.getAffineCover();
	}

	@Override
	public List<Integer> affineCoverIndex(ProjectivePoint<T> p) {
		return asProjectiveVariety.affineCoverIndex(p);
	}

	@Override
	public AffinePoint<T> asAffinePoint(ProjectivePoint<T> p, int affineCoverIndex) {
		return asProjectiveVariety.asAffinePoint(p, affineCoverIndex);
	}

	@Override
	public ProjectivePoint<T> fromAffinePoint(AffinePoint<T> p, int affineCoverIndex) {
		return asProjectiveVariety.fromAffinePoint(p, affineCoverIndex);
	}

	@Override
	public Morphism<T, AffinePoint<T>, ProjectivePoint<T>> embedding(int affineCoverIndex) {
		return asProjectiveVariety.embedding(affineCoverIndex);
	}

	@Override
	public Morphism<T, SingletonPoint, ProjectivePoint<T>> pointAsMorphism(ProjectivePoint<T> p) {
		return asProjectiveVariety.pointAsMorphism(p);
	}

}
