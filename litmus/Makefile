#
# Makefile for LITMUS^RT
#

obj-y     = sched_plugin.o litmus.o \
	    preempt.o \
	    litmus_proc.o \
	    budget.o \
	    clustered.o \
	    jobs.o \
	    sync.o \
	    rt_domain.o \
	    edf_common.o \
	    fp_common.o \
	    fdso.o \
	    locking.o \
	    srp.o \
	    bheap.o \
	    binheap.o \
	    ctrldev.o \
	    uncachedev.o \
	    sched_gsn_edf.o \
	    sched_psn_edf.o \
#	    sched_pfp.o

obj-$(CONFIG_PLUGIN_CEDF) += sched_cedf.o
obj-$(CONFIG_PLUGIN_PFAIR) += sched_pfair.o

obj-$(CONFIG_FEATHER_TRACE) += ft_event.o ftdev.o
obj-$(CONFIG_SCHED_TASK_TRACE) += sched_task_trace.o
obj-$(CONFIG_SCHED_DEBUG_TRACE) += sched_trace.o
obj-$(CONFIG_SCHED_OVERHEAD_TRACE) += trace.o

obj-y += sched_pres.o

obj-y += reservations/
obj-y += sched_demo.o
obj-y += simple_heap.o
