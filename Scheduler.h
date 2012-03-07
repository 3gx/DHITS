#ifndef _SCHEDULER_
#define _SCHEDULER_

template<int RUNGMAX>
struct Scheduler 
{
//	enum {RUNGMAX = 32};

	double dtmax;
	double dt_tick;
	unsigned long long tsysU;
	int min_rung, max_rung;     
	std::vector<int> list[RUNGMAX];

	Scheduler() {}
	~Scheduler() {};

	Scheduler(const double _dtmax) : dtmax(_dtmax) 
	{
		for (int i = 0; i < RUNGMAX; i++) {
			list[i].reserve(64);
			list[i].clear();
		}
		dt_tick = dtmax / (1LLU << (RUNGMAX - 1));
		fprintf(stderr, " dt_tick= %lg\n", dt_tick);
		tsysU = 0;
		min_rung = 0;
		max_rung = RUNGMAX - 1;
	}

	// IEEE complient ONLY !!!!!!!!!
	static const int exp_of(const double dt) 
	{
		union {double d; unsigned long long u;} pack;
		pack.d = dt;
		return (pack.u >> 52) & 0xfff;
	};

	void flush_list() 
  {
		for (int i = 0; i < RUNGMAX; i++)
      list[i].clear();
	}

	const int push_particle(const int index, const int rung) 
	{
		assert(rung < RUNGMAX);
		list[rung].push_back(index);
		return rung;
	}

	const int push_particle(const int index, const double dt) 
	{
		const int rung = std::max(exp_of(dtmax) - exp_of(dt), min_rung);
		assert(rung < RUNGMAX);
		list[rung].push_back(index);
		return rung; //dtmax/(1U << rung);
	};


	const int get_rung(const double dt) 
	{
		const int rung = std::max(exp_of(dtmax) - exp_of(dt), min_rung);
		assert(rung < RUNGMAX);
		return rung; 
	}

	const int get_rung(const double dt, const int _rung) 
	{
		const int rung = std::max(exp_of(dtmax) - exp_of(dt), _rung);
		assert(rung < RUNGMAX);
		return rung; 
	}

	const double get_dt(const int rung) {
		unsigned long long dtU = 1LLU << (RUNGMAX - 1 - rung);
		return dt_tick * dtU;
	}

	const double pull_active_list(std::vector<int> &ptcl_list)
	{
		max_rung= 0;
		for (int i = 0; i < RUNGMAX; i++) 
			if (list[i].size() > 0) max_rung= i;

		assert(max_rung != 31);

		const unsigned long long dtU    = 1LLU << (RUNGMAX - 1 - max_rung);
		const unsigned long long tnextU = (tsysU/dtU + 1)*dtU;
		const unsigned long long flags  = tsysU ^ tnextU;

		for (int i = 0; i < RUNGMAX; i++) {
			if (flags & (1LLU << i)) {
				if (list[RUNGMAX - 1 - i].size() > 0) { 
					for (size_t j = 0; j < (const size_t)list[RUNGMAX - 1 - i].size(); j++) {
						const int idx = list[RUNGMAX - 1 - i][j];
						ptcl_list.push_back(idx);
					}
					list[RUNGMAX - 1 - i].clear();
				}
				min_rung = RUNGMAX - 1 - i;
			}
		}
    const unsigned long long tprevU = tsysU;
		tsysU = tnextU;
		
    return dt_tick * (tsysU - tprevU);
	}

	void debug_dump() const 
  { 
		FILE *ferr = stderr;
		int n_tot = 0;
		for (int i = 0; i < RUNGMAX; i++) {
			const int count = list[i].size();
			n_tot += count;
			fprintf(ferr, " %4d : %6d \n", i, count);
		}
		fprintf(ferr, " total: %6d \n",  n_tot);
	}

	const double get_tsys() const 
  {
		return dt_tick * tsysU;
	}
  void set_tsys(const double t) 
  {
    assert(fmod(t, dt_tick) == 0.0);
    tsysU = (unsigned long long)(t/dt_tick);
  }


  int get_n() const 
  {
    int n = 0;
    for (int i = 0; i < RUNGMAX; i++)
      n += list[i].size();
    return n;
  }

	const double get_dt_pred(const int rung) const 
  {
		const unsigned long long dtU = 1LLU << (RUNGMAX - 1 - rung);
		const unsigned long long tmp = (rung >= min_rung) ? dtU : tsysU & (dtU - 1);
		return dt_tick * tmp;
	}

	const double get_dt_corr(const int rung) const 
  {
		return dtmax / (1LLU << rung) ;
	}

};

#endif
