#include <linux/sched.h>
#include <linux/percpu.h>
#include <linux/hrtimer.h>
#include <linux/uaccess.h>
#include <linux/module.h>

#include <litmus/debug_trace.h>
#include <litmus/litmus.h>
#include <litmus/preempt.h>
#include <litmus/sched_plugin.h>
#include <litmus/np.h>

#include <litmus/budget.h>

//--SS--
#define ns2ms(x) ((x) / 1000000)

struct enforcement_timer {
	/* The enforcement timer is used to accurately police
	 * slice budgets. */
	struct hrtimer		timer;
	int			armed;
};

DEFINE_PER_CPU(struct enforcement_timer, budget_timer);

static enum hrtimer_restart on_enforcement_timeout(struct hrtimer *timer)
{
	struct enforcement_timer* et = container_of(timer,
						    struct enforcement_timer,
						    timer);
	unsigned long flags;

	local_irq_save(flags);
	TRACE("enforcement timer for LC fired.\n");
	et->armed = 0;
	/* activate scheduler */
	litmus_reschedule_local();
	local_irq_restore(flags);

	return  HRTIMER_NORESTART;
}

// --SS--
static enum hrtimer_restart on_enforcement_timeout_for_HC(struct hrtimer *timer)
{
	struct enforcement_timer* et = container_of(timer,
						    struct enforcement_timer,
						    timer);
	unsigned long flags;

	local_irq_save(flags);
  /*TRACE("Fired at %llu. HC Timer was added at %llu for pid %d\n",
  litmus_clock(), et->timer.when_added,  et->timer.start_pid);*/
  // --SS-- increment the System mode here
  //TRACE("System Mode was: %u\n", get_system_mode());
  if(!IS_LO_SYSTEM_MODE) {
    TRACE("Enforcement timer for HC fired in System Mode 1 --SS--\n");
  } else {
    TRACE("Enforcement timer for HC fired in System Mode 0:--SS--\n");
    change_system_mode(1);
  }
  //TRACE("System Mode is now: %u\n", get_system_mode());
	et->armed = 0;
	/* activate scheduler */
	//litmus_reschedule_local();
	local_irq_restore(flags);

	return  HRTIMER_NORESTART;
}

/* assumes called with IRQs off */
static void cancel_enforcement_timer(struct enforcement_timer* et)
{
	int ret;

	TRACE("cancelling enforcement timer.\n");

	/* Since interrupts are disabled and et->armed is only
	 * modified locally, we do not need any locks.
	 */

	if (et->armed) {
		ret = hrtimer_try_to_cancel(&et->timer);
		/* Should never be inactive. */
		BUG_ON(ret == 0);
		/* Should never be running concurrently. */
		BUG_ON(ret == -1);

		et->armed = 0;
	}
}

/* --SS-- this timer is for HC task's monitoring of LO-time.
This is used to change the System_Mode. This function is only to enable the
timer but not the actual timer itself.
*/
static void arm_enforcement_timer_for_HC(struct enforcement_timer* et,
  struct task_struct* t) {
  lt_t when_to_fire;

	WARN_ONCE(!hrtimer_is_hres_active(&et->timer),
		KERN_ERR "WARNING: no high resolution timers available for HC!?\n");

	/* Calling this when there is no budget left for the task
	 * makes no sense, unless the task is non-preemptive. */
	BUG_ON(budget_exhausted(t) && (!is_np(t)));

	/* hrtimer_start_range_ns() cancels the timer
	 * anyway, so we don't have to check whether it is still armed */

	if (likely(!is_np(t))) {
		when_to_fire = litmus_clock() + budget_remaining(t);
		et->timer.function = on_enforcement_timeout_for_HC;
    /*TRACE_TASK(t, "arming enforcement timer for HC, budget_rem: %llu at %llu;"
    "will fire at %llu\n",
    budget_remaining(t), litmus_clock(), when_to_fire);*/
    /*et->timer.start_pid = (t)->pid;
    et->timer.when_added = ns_to_ktime(clock_current);*/
		hrtimer_start(&et->timer, ns_to_ktime(when_to_fire),
			HRTIMER_MODE_ABS_PINNED);
		et->armed = 1;
	}
  //TRACE_TASK(t, "timer duration: %llu\n", when_to_fire - litmus_clock());
}

void elongate_enforcement_timer(struct task_struct* t, lt_t extra_time) {
  lt_t when_to_fire;
	struct enforcement_timer* et = this_cpu_ptr(&budget_timer);
  if(!t) {
    TRACE_TASK(t, "NULL task to elongate\n");
  }
  cancel_enforcement_timer(et);
  when_to_fire = litmus_clock() + budget_remaining(t) + extra_time;
  et->timer.function = on_enforcement_timeout_for_HC;
  hrtimer_start(&et->timer, ns_to_ktime(when_to_fire),
    HRTIMER_MODE_ABS_PINNED);
  et->armed = 1;
  TRACE_TASK(t, "elongated HC timer by %llu\n", ns2ms(extra_time));
}

static void local_reschedule_call_for_LC_timeout(struct enforcement_timer* et) {
  // taken from the above enforcement time out function
	unsigned long flags;
	local_irq_save(flags);
	TRACE("enforcement timer for LC fired.\n");
  cancel_enforcement_timer(et);
	/* activate scheduler */
	litmus_reschedule_local();
	local_irq_restore(flags);
}

/* assumes called with IRQs off */
static void arm_enforcement_timer(struct enforcement_timer* et,
				  struct task_struct* t)
{
	lt_t when_to_fire;
  TRACE_TASK(t, "arming enforcement timer for LC at %llu\n",
  litmus_clock());

	WARN_ONCE(!hrtimer_is_hres_active(&et->timer),
		KERN_ERR "WARNING: no high resolution timers available!?\n");

  TRACE_TASK(t, "LC task budget rem: %llu, Mode: %u\n", budget_remaining(t),
  get_system_mode());
	if(budget_exhausted(t) && (!is_np(t))) {
    TRACE_TASK(t, "Calling timer timeout for LC budget exhaustion\n");
    local_reschedule_call_for_LC_timeout(et);
    return;
  }
	
  /* Calling this when there is no budget left for the task
	 * makes no sense, unless the task is non-preemptive. */
  BUG_ON(budget_exhausted(t) && (!is_np(t)));

	/* hrtimer_start_range_ns() cancels the timer
	 * anyway, so we don't have to check whether it is still armed */

	if (likely(!is_np(t))) {
		when_to_fire = litmus_clock() + budget_remaining(t);
		et->timer.function = on_enforcement_timeout;
		hrtimer_start(&et->timer, ns_to_ktime(when_to_fire),
			HRTIMER_MODE_ABS_PINNED);
		et->armed = 1;
	}
}


/* expects to be called with IRQs off */
void update_enforcement_timer(struct task_struct* t)
{
	struct enforcement_timer* et = this_cpu_ptr(&budget_timer);

	if (t && budget_precisely_enforced(t)) {
		/* Make sure we call into the scheduler when this budget
		 * expires. */
    // --SS-- call different timer enforcement functions for different class of
    // task
    if(is_hrt(t))
		  arm_enforcement_timer_for_HC(et, t);
    else
      arm_enforcement_timer(et, t);
	} else if (et->armed) {
		/* Make sure we don't cause unnecessary interrupts. */
		cancel_enforcement_timer(et);
	}
}


static int __init init_budget_enforcement(void)
{
	int cpu;
	struct enforcement_timer* et;

	for (cpu = 0; cpu < NR_CPUS; cpu++)  {
		et = &per_cpu(budget_timer, cpu);
		hrtimer_init(&et->timer, CLOCK_MONOTONIC, HRTIMER_MODE_ABS);
		et->timer.function = on_enforcement_timeout;
	}
	return 0;
}

void litmus_current_budget(lt_t *used_so_far, lt_t *remaining)
{
	struct task_struct *t = current;
	unsigned long flags;
	s64 delta;

	local_irq_save(flags);

	delta = sched_clock_cpu(smp_processor_id()) - t->se.exec_start;
	if (delta < 0)
		delta = 0;

	/*TRACE_CUR("current_budget: sc:%llu start:%llu lt_t:%llu delta:%lld exec-time:%llu rem:%llu\n",
		sched_clock_cpu(smp_processor_id()), t->se.exec_start,
		litmus_clock(), delta,
		tsk_rt(t)->job_params.exec_time,
		budget_remaining(t));*/

	if (used_so_far)
		*used_so_far = tsk_rt(t)->job_params.exec_time + delta;

	if (remaining) {
		*remaining = budget_remaining(t);
		if (*remaining > delta)
			*remaining -= delta;
		else
			*remaining = 0;
	}

	local_irq_restore(flags);
}

asmlinkage long sys_get_current_budget(
	lt_t __user * _expended,
	lt_t __user *_remaining)
{
	lt_t expended = 0, remaining = 0;

	if (is_realtime(current))
		litmus->current_budget(&expended, &remaining);

	if (_expended && put_user(expended, _expended))
		return -EFAULT;

	if (_remaining && put_user(remaining, _remaining))
		return -EFAULT;

	return 0;
}

module_init(init_budget_enforcement);
