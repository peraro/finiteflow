#include <fflow/mp_gcd.hh>
#include <fflow/gcd.hh>

namespace fflow {

  bool rat_mod(const MPRational & q, const MPInt & p,
               MPInt & res)
  {
    MPInt iden;
    bool ret = mpz_invert(iden.get(), mpq_denref(q.get()), p.get());
    if (ret) {
      mpz_mul(res.get(), mpq_numref(q.get()), iden.get());
      mpz_mod(res.get(), res.get(), p.get());
      return true;
    } else {
      return false;
    }
  }

  bool rat_rec(const MPInt & a, const MPInt & b,
               MPRational & res)
  {
    //Int t = 1, tt = 0;
    MPInt s, ss(1), r(b), rr(a);
    MPInt q, tmp;

    // rr2 = rr^2
    MPInt rr2;
    mpz_mul(rr2.get(), rr.get(), rr.get());

    while (rr2.cmp(b) > 0) {

      if (r.sign() == 0)
        return false;

      // q = rr/r
      mpz_div(q.get(), rr.get(), r.get());

      // tmp = rr;  rr = r;  r = tmp - q * r;
      tmp = std::move(rr);
      mpz_submul(tmp.get(), q.get(), r.get());
      rr = std::move(tmp);
      rr.swap(r);

      // tmp = ss;  ss = s;  s = tmp - q * s;
      tmp = std::move(ss);
      mpz_submul(tmp.get(), q.get(), s.get());
      ss = std::move(tmp);
      ss.swap(s);

      // rr2 == rr^2
      mpz_mul(rr2.get(), rr.get(), rr.get());
    }

    mpz_swap(mpq_numref(res.get()), rr.get());
    mpz_swap(mpq_denref(res.get()), ss.get());

    res.canonicalize();

    return true;
  }

  void chinese_remainder_coeffs(const MPInt & p, UInt q,
                                MPInt & c1, MPInt & c2, MPInt & mod12)
  {
    MPInt qp, pinvp, qinvq;
    Mod modq(q);
    UInt ipmodq;

    // qp = q*p;
    mpz_mul_ui(qp.get(), p.get(), q);

    // ipmodq = mul_inv(p % q, q)
    ipmodq = mpz_fdiv_ui(p.get(), q);
    ipmodq = mul_inv(ipmodq, modq);

    // pinvp = p * ipmodq % qp
    mpz_mul_ui(pinvp.get(), p.get(), ipmodq);
    mpz_mod(pinvp.get(), pinvp.get(), qp.get());

    // qinvq = (1-pinvp) % qp
    mpz_ui_sub(qinvq.get(), 1, pinvp.get());
    mpz_mod(qinvq.get(), qinvq.get(), qp.get());

    // store the results
    qp.swap(mod12);
    qinvq.swap(c1);
    pinvp.swap(c2);
  }

  void chinese_remainder_from_coeffs(const MPInt & a, UInt b,
                                     const MPInt & c1, const MPInt & c2,
                                     const MPInt & mod12,
                                     MPInt & res)
  {
    MPInt r;
    mpz_mul_ui(r.get(), c2.get(), b);
    mpz_addmul(r.get(), c1.get(), a.get());
    mpz_mod(r.get(), r.get(), mod12.get());
    r.swap(res);
  }

} // namespace fflow
