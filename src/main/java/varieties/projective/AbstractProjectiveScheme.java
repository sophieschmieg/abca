package varieties.projective;

import java.util.List;

import fields.interfaces.Element;
import varieties.AbstractScheme;
import varieties.FunctionField;
import varieties.Morphism;
import varieties.SpectrumOfField.SingletonPoint;
import varieties.affine.AffineCover;
import varieties.affine.AffinePoint;

public abstract class AbstractProjectiveScheme<T extends Element<T>> extends AbstractScheme<T, ProjectivePoint<T>>
		implements ProjectiveScheme<T> {
	private GenericProjectiveScheme<T> asGenericProjectiveScheme;
	private FunctionField<T> functionField;

	public AbstractProjectiveScheme(GenericProjectiveScheme<T> asProjectiveVariety) {
		this.asGenericProjectiveScheme = asProjectiveVariety;
	}

	@Override
	public Exactness exactness() {
		return asGenericProjectiveScheme.exactness();
	}

	@Override
	public GenericProjectiveScheme<T> asGenericProjectiveScheme() {
		return asGenericProjectiveScheme;
	}

	@Override
	public AffineCover<T> getAffineCover() {
		return asGenericProjectiveScheme.getAffineCover();
	}

	@Override
	public List<Integer> affineCoverIndex(ProjectivePoint<T> p) {
		return asGenericProjectiveScheme.affineCoverIndex(p);
	}

	@Override
	public AffinePoint<T> asAffinePoint(ProjectivePoint<T> p, int affineCoverIndex) {
		return asGenericProjectiveScheme.asAffinePoint(p, affineCoverIndex);
	}

	@Override
	public ProjectivePoint<T> fromAffinePoint(AffinePoint<T> p, int affineCoverIndex) {
		return asGenericProjectiveScheme.fromAffinePoint(p, affineCoverIndex);
	}

	@Override
	public Morphism<T, AffinePoint<T>, ProjectivePoint<T>> embedding(int affineCoverIndex) {
		return asGenericProjectiveScheme.embedding(affineCoverIndex);
	}

	@Override
	public Morphism<T, SingletonPoint, ProjectivePoint<T>> pointAsMorphism(ProjectivePoint<T> p) {
		return asGenericProjectiveScheme.pointAsMorphism(p);
	}

	@Override
	public final FunctionField<T> getFunctionField() {
		if (functionField == null) {
			functionField = new FunctionField<>(this);
		}
		return functionField;
	}

}
