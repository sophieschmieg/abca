package varieties;

import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.MathSet;
import varieties.SpectrumOfField.SingletonPoint;
import varieties.affine.AffineCover;
import varieties.affine.AffinePoint;

public interface Scheme<T extends Element<T>, S extends Element<S>> extends MathSet<S> {
	public Field<T> getField();
	public boolean hasRationalPoint(S p);
	public int dimension();
	public AffineCover<T> getAffineCover();
	public List<Integer> affineCoverIndex(S p);
	public AffinePoint<T> asAffinePoint(S p, int affineCoverIndex);
	public S fromAffinePoint(AffinePoint<T> p, int affineCoverIndex);
	public Morphism<T, AffinePoint<T>, S> embedding(int affineCoverIndex);
	public Morphism<T, SingletonPoint, S> pointAsMorphism(S p);
	public boolean isConnected();
	public boolean isIntegral();
	public boolean isReduced();
	public boolean isIrreducible();
	public List<? extends Scheme<T, S>> irreducibleComponents();
}
