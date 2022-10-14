#ifndef _FTK_FEATURE_CURVE_SET_PP_HH
#define _FTK_FEATURE_CURVE_SET_PP_HH

#include <ftk/config.hh>
#include <ftk/features/feature_curve_set.hh>
#include <ftk/utils/string.hh>

namespace ftk {

struct feature_curve_set_post_processor_t {
  feature_curve_set_post_processor_t(const std::string str) { parse_ops(str); }

  void parse_ops(const std::string str) { ops = split(str, ","); }
  void reset() { ops.clear(); }

  void filter(feature_curve_set_t& trajs) const;
 
protected:
  std::vector<std::string> ops;
};

inline void feature_curve_set_post_processor_t::filter(feature_curve_set_t &trajs) const
{
  for (const auto op : ops) {
    // if (op == "discard_high_cond") {
    //   trajs.foreach([](ftk::feature_curve_t& t) {
    //     t.discard_high_cond();
    //     t.update_statistics();
    //   });
    // } else 
    if (op == "smooth_types") {
      trajs.foreach([](ftk::feature_curve_t& t) {
        t.smooth_ordinal_types();
        t.smooth_interval_types();
        t.update_statistics();
      });
    } else if (op == "rotate") {
      trajs.foreach([](ftk::feature_curve_t& t) {
        t.rotate();
        t.update_statistics();
      });
    } else if (starts_with(op, "duration_pruning")) {
      // syntax: duration_pruning:value, e,g. duration_pruning:2.0
      // TODO: parse duration threshold
    } else if (op == "split") {
      trajs.split_all();
    } else if (op == "discard_interval_points") {
      trajs.foreach([](ftk::feature_curve_t& t) {
        t.discard_interval_points();
        t.update_statistics();
      });
    } else if (op == "reorder") {
      trajs.foreach([](ftk::feature_curve_t& t) {
        t.reorder();
        t.update_statistics();
      });
    } else if (op == "adjust_time") {
      trajs.foreach([](ftk::feature_curve_t& t) {
        t.adjust_time();
        t.update_statistics();
      });
    } else if (op == "derive_velocity") {
      trajs.foreach([](ftk::feature_curve_t& t) {
        t.discard_interval_points();
        t.derive_velocity();
        t.update_statistics();
      });
    } else 
      fatal(FTK_ERR_UNKNOWN_OPTIONS);
  }
}

} // namespace ftk

#endif
