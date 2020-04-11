#ifndef _LITMUS_BUDGET_H_
#define _LITMUS_BUDGET_H_

// --SS-- if exec_cost is zero -- this happens to LC task in HI-mode
#define is_zero_exec_cost(t) (get_exec_cost(t) == 0)

/* Update the per-processor enforcement timer (arm/reproram/cancel) for
 * the next task. */
void update_enforcement_timer(struct task_struct* t);

// The following elongates an enforcement timer of a task by a specified time
void elongate_enforcement_timer(struct task_struct* t, lt_t extra_time);

inline static int budget_exhausted(struct task_struct* t)
{
  // zero exec cost happens to LC task when the System Mode is HI
  if(is_zero_exec_cost(t))
    return 1;
  if(IS_LO_SYSTEM_MODE && is_elongated(t))
    return get_exec_time(t) >= get_elongated_exec_cost(t);
	return get_exec_time(t) >= get_exec_cost(t);
}

inline static lt_t budget_remaining(struct task_struct* t)
{
	if (!budget_exhausted(t)) {
    if(IS_LO_SYSTEM_MODE && is_elongated(t)) {
      return get_elongated_exec_cost(t) - get_exec_time(t);
    }
		return get_exec_cost(t) - get_exec_time(t);
  }
	else
		/* avoid overflow */
		return 0;
}

#define budget_enforced(t) (tsk_rt(t)->task_params.budget_policy != NO_ENFORCEMENT)

#define budget_precisely_enforced(t) (tsk_rt(t)->task_params.budget_policy \
				      == PRECISE_ENFORCEMENT)

static inline int requeue_preempted_job(struct task_struct* t)
{
	/* Add task to ready queue only if not subject to budget enforcement or
	 * if the job has budget remaining. t may be NULL.
	 */
	return t && !is_completed(t) &&
		(!budget_exhausted(t) || !budget_enforced(t));
}

void litmus_current_budget(lt_t *used_so_far, lt_t *remaining);

#endif
