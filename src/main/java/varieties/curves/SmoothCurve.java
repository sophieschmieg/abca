package varieties.curves;

import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.polynomials.CoordinateRing;
import varieties.RationalFunction;
import varieties.Scheme;
import varieties.curves.DivisorGroup.Divisor;
import varieties.projective.ProjectivePoint;
import varieties.projective.GenericProjectiveScheme;

public interface SmoothCurve<T extends Element<T>> extends Scheme<T, ProjectivePoint<T>> {
	public int getEmbeddingDimension();
	public GenericProjectiveScheme<T> asGenericProjectiveScheme();
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p);
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p);
	public List<RationalFunction<T>> getRiemannRochSpace(Divisor<T> div);
	public boolean isPrincipal(Divisor<T> div);
	public CoordinateRing<T> getCoordinateRing();
}
