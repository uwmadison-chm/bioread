# GENERATED FILE -- do not edit.
# Regenerate with: python notes/kaitai/gen_ksy.py > notes/kaitai/acq.ksy
#
# A Kaitai Struct description of the AcqKnowledge header layouts, generated
# from the _versioned_fields tables in bioread/headers.py.
#
# This describes the headers only. Locating them in the file -- summing channel
# data sizes, scanning for the channel dtype headers -- is not expressible in
# Kaitai's expression language and stays in Python. See README.md.

meta:
  id: acq
  file-extension: acq

instances:
  # bioread works out the byte order by reading the version field both
  # ways and keeping whichever gives the smaller number. That heuristic
  # is expressible here, so the spec is self-orienting.
  version_le:
    pos: 2
    type: u4le
  version_be:
    pos: 2
    type: u4be
  is_le:
    value: version_le <= version_be
  file_revision:
    value: 'is_le ? version_le : version_be'
  graph_header:
    pos: 0
    type:
      switch-on: file_revision >= 61
      cases:
        'true': graph_header_post4
        'false': graph_header_pre4

types:
  graph_header_pre4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: n_item_header_len
        type: s2
      - id: l_version
        type: u4
      - id: l_ext_item_header_len
        type: s4
        if: _root.file_revision >= 30
      - id: n_channels
        type: s2
        if: _root.file_revision >= 30
      - id: n_horiz_axis_type
        type: s2
        if: _root.file_revision >= 30
      - id: n_cur_channel
        type: s2
        if: _root.file_revision >= 30
      - id: d_sample_time
        type: f8
        if: _root.file_revision >= 30
      - id: d_time_offset
        type: f8
        if: _root.file_revision >= 30
      - id: d_time_scale
        type: f8
        if: _root.file_revision >= 30
      - id: d_time_cursor1
        type: f8
        if: _root.file_revision >= 30
      - id: d_time_cursor2
        type: f8
        if: _root.file_revision >= 30
      - id: rc_window
        type: s2
        repeat: expr
        repeat-expr: 4
        if: _root.file_revision >= 30
      - id: n_measurement
        type: s2
        repeat: expr
        repeat-expr: 6
        if: _root.file_revision >= 30
      - id: f_hilite
        type: s2
        if: _root.file_revision >= 30
      - id: d_first_time_offset
        type: f8
        if: _root.file_revision >= 30
      - id: n_rescale
        type: s2
        if: _root.file_revision >= 30
      - id: sz_horiz_units1
        size: 40
        if: _root.file_revision >= 30
      - id: sz_horiz_units2
        size: 10
        if: _root.file_revision >= 30
      - id: n_in_memory
        type: s2
        if: _root.file_revision >= 30
      - id: f_grid
        type: s2
        if: _root.file_revision >= 30
      - id: f_markers
        type: s2
        if: _root.file_revision >= 30
      - id: n_plot_draft
        type: s2
        if: _root.file_revision >= 30
      - id: n_disp_mode
        type: s2
        if: _root.file_revision >= 30
      - id: r_r_reserved
        type: s2
        if: _root.file_revision >= 30
      - id: b_show_tool_bar
        type: s2
        if: _root.file_revision >= 34
      - id: b_show_channel_buttons
        type: s2
        if: _root.file_revision >= 34
      - id: b_show_measurements
        type: s2
        if: _root.file_revision >= 34
      - id: b_show_markers
        type: s2
        if: _root.file_revision >= 34
      - id: b_show_journal
        type: s2
        if: _root.file_revision >= 34
      - id: cur_x_channel
        type: s2
        if: _root.file_revision >= 34
      - id: mmt_precision
        type: s2
        if: _root.file_revision >= 34
      - id: n_measurement_rows
        type: s2
        if: _root.file_revision >= 35
      - id: mmt40
        type: s2
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 35
      - id: mmt_chan40
        type: s2
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 35
      - id: mmt_calc_opnd1
        type: s2
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 36
      - id: mmt_calc_opnd2
        type: s2
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 36
      - id: mmt_calc_op
        type: s2
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 36
      - id: mmt_calc_constant
        type: f8
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 36
      - id: b_new_grid_with_minor
        type: s4
        if: _root.file_revision >= 38
      - id: color_major_grid
        size: 4
        if: _root.file_revision >= 38
      - id: color_minor_grid
        size: 4
        if: _root.file_revision >= 38
      - id: w_major_grid_style
        type: u2
        if: _root.file_revision >= 38
      - id: w_minor_grid_style
        type: u2
        if: _root.file_revision >= 38
      - id: w_major_grid_width
        type: u2
        if: _root.file_revision >= 38
      - id: w_minor_grid_width
        type: u2
        if: _root.file_revision >= 38
      - id: b_fixed_units_div
        type: u4
        if: _root.file_revision >= 38
      - id: b_mid__range__show
        type: u4
        if: _root.file_revision >= 38
      - id: d_start__middle__point
        type: f8
        if: _root.file_revision >= 38
      - id: d_offset__point
        type: f8
        repeat: expr
        repeat-expr: 60
        if: _root.file_revision >= 38
      - id: h_grid
        type: f8
        if: _root.file_revision >= 38
      - id: v_grid
        type: f8
        repeat: expr
        repeat-expr: 60
        if: _root.file_revision >= 38
      - id: b_enable_wave_tools
        type: s4
        if: _root.file_revision >= 38
      - id: hoziz_precision
        type: s2
        if: _root.file_revision >= 39
      - id: reserved
        size: 20
        if: _root.file_revision >= 41
      - id: b_overlap_mode
        type: s4
        if: _root.file_revision >= 41
      - id: b_show_hardware
        type: s4
        if: _root.file_revision >= 41
      - id: b_x_auto_plot
        type: s4
        if: _root.file_revision >= 41
      - id: b_x_auto_scroll
        type: s4
        if: _root.file_revision >= 41
      - id: b_start_button_visible
        type: s4
        if: _root.file_revision >= 41
      - id: b_compressed
        type: s4
        if: _root.file_revision >= 41
      - id: b_always_start_button_visible
        type: s4
        if: _root.file_revision >= 41
      - id: path_video
        size: 260
        if: _root.file_revision >= 43
      - id: opt_sync_delay
        type: s4
        if: _root.file_revision >= 43
      - id: sync_delay
        type: f8
        if: _root.file_revision >= 43
      - id: b_hrp__paste_measurements
        type: s4
        if: _root.file_revision >= 43
      - id: graph_type
        type: s4
        if: _root.file_revision >= 45
      - id: mmt_calc_expr
        size: 10240
        if: _root.file_revision >= 45
      - id: mmt_moment_order
        type: s4
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 45
      - id: mmt_time_delay
        type: s4
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 45
      - id: mmt_embed_dim
        type: s4
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 45
      - id: mmt_mi_delay
        type: s4
        repeat: expr
        repeat-expr: 40
        if: _root.file_revision >= 45

  graph_header_post4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: n_item_header_len
        type: s2
      - id: l_version
        type: s4
      - id: l_ext_item_header_len
        type: s4
        if: _root.file_revision >= 30
      - id: n_channels
        type: s2
        if: _root.file_revision >= 30
      - id: n_horiz_axis_type
        type: s2
        if: _root.file_revision >= 30
      - id: n_cur_channel
        type: s2
        if: _root.file_revision >= 30
      - id: d_sample_time
        type: f8
        if: _root.file_revision >= 30
      - id: d_time_offset
        type: f8
        if: _root.file_revision >= 30
      - id: d_time_scale
        type: f8
        if: _root.file_revision >= 30
      - id: d_time_cursor1
        type: f8
        if: _root.file_revision >= 30
      - id: d_time_cursor2
        type: f8
        if: _root.file_revision >= 30
      - id: rc_window
        type: s2
        repeat: expr
        repeat-expr: 4
        if: _root.file_revision >= 30
      - id: n_measurement
        type: s2
        repeat: expr
        repeat-expr: 6
        if: _root.file_revision >= 30
      - id: f_hilite
        type: s2
        if: _root.file_revision >= 30
      - id: d_first_time_offset
        type: f8
        if: _root.file_revision >= 30
      - id: n_rescale
        type: s2
        if: _root.file_revision >= 30
      - id: sz_horiz_units1
        size: 40
        if: _root.file_revision >= 30
      - id: sz_horiz_units2
        size: 10
        if: _root.file_revision >= 30
      - id: n_in_memory
        type: s2
        if: _root.file_revision >= 30
      - id: f_grid
        type: s2
        if: _root.file_revision >= 30
      - id: f_markers
        type: s2
        if: _root.file_revision >= 30
      - id: n_plot_draft
        type: s2
        if: _root.file_revision >= 30
      - id: n_disp_mode
        type: s2
        if: _root.file_revision >= 30
      - id: r_r_reserved
        type: s2
        if: _root.file_revision >= 30
      - id: unknown
        size: 822
        if: _root.file_revision >= 61
      - id: b_compressed
        type: s4
        if: _root.file_revision >= 61
      - id: unknown2
        size: 1422
        if: _root.file_revision >= 61
      - id: h_expected_paddings
        type: s2
        if: _root.file_revision >= 124

  unknown_padding_header:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_channel_len
        type: s4
      - id: uknown
        size: 36

  channel_header_pre4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_chan_header_len
        type: s4
        if: _root.file_revision >= 30
      - id: n_num
        type: s2
        if: _root.file_revision >= 30
      - id: sz_comment_text
        size: 40
        if: _root.file_revision >= 30
      - id: rgb_color
        size: 4
        if: _root.file_revision >= 30
      - id: n_disp_chan
        type: s2
        if: _root.file_revision >= 30
      - id: d_volt_offset
        type: f8
        if: _root.file_revision >= 30
      - id: d_volt_scale
        type: f8
        if: _root.file_revision >= 30
      - id: sz_units_text
        size: 20
        if: _root.file_revision >= 30
      - id: l_buf_length
        type: s4
        if: _root.file_revision >= 30
      - id: d_ampl_scale
        type: f8
        if: _root.file_revision >= 30
      - id: d_ampl_offset
        type: f8
        if: _root.file_revision >= 30
      - id: n_chan_order
        type: s2
        if: _root.file_revision >= 30
      - id: n_disp_size
        type: s2
        if: _root.file_revision >= 30
      - id: plot_mode
        type: s2
        if: _root.file_revision >= 34
      - id: v_mid
        type: f8
        if: _root.file_revision >= 34
      - id: sz_description
        size: 128
        if: _root.file_revision >= 38
      - id: n_var_sample_divider
        type: s2
        if: _root.file_revision >= 38
      - id: vert_precision
        type: s2
        if: _root.file_revision >= 39
      - id: active_segment_color
        size: 4
        if: _root.file_revision >= 43
      - id: active_segment_style
        type: s4
        if: _root.file_revision >= 43

  channel_header_post4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_chan_header_len
        type: s4
        if: _root.file_revision >= 30
      - id: n_num
        type: s2
        if: _root.file_revision >= 30
      - id: sz_comment_text
        size: 40
        if: _root.file_revision >= 30
      - id: not_color
        size: 4
        if: _root.file_revision >= 30
      - id: n_disp_chan
        type: s2
        if: _root.file_revision >= 30
      - id: d_volt_offset
        type: f8
        if: _root.file_revision >= 30
      - id: d_volt_scale
        type: f8
        if: _root.file_revision >= 30
      - id: sz_units_text
        size: 20
        if: _root.file_revision >= 30
      - id: l_buf_length
        type: s4
        if: _root.file_revision >= 30
      - id: d_ampl_scale
        type: f8
        if: _root.file_revision >= 30
      - id: d_ampl_offset
        type: f8
        if: _root.file_revision >= 30
      - id: n_chan_order
        type: s2
        if: _root.file_revision >= 30
      - id: n_disp_size
        type: s2
        if: _root.file_revision >= 30
      - id: unknown
        size: 40
        if: _root.file_revision >= 61
      - id: n_var_sample_divider
        type: s2
        if: _root.file_revision >= 61

  foreign_header_pre4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: n_length
        type: s2
        if: _root.file_revision >= 30
      - id: n_type
        type: s2
        if: _root.file_revision >= 30

  foreign_header_post4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_length
        type: s4
        if: _root.file_revision >= 61

  channel_d_type_header:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: n_size
        type: s2
        if: _root.file_revision >= 30
      - id: n_type
        type: s2
        if: _root.file_revision >= 30

  journal_header_pre4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: tag
        size: 4
        if: _root.file_revision >= 30
      - id: h_show
        type: s2
        if: _root.file_revision >= 30
      - id: l_journal_len
        type: s4
        if: _root.file_revision >= 30

  journal_header_post4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: b_unknown1
        size: 262
        if: _root.file_revision >= 61
      - id: l_early_journal_len
        type: s4
        if: _root.file_revision >= 61
      - id: b_unknown2
        size: 290
        if: _root.file_revision >= 61
      - id: b_unknown3
        size: 26
        if: _root.file_revision >= 108
      - id: b_unknown4
        size: 4
        if: _root.file_revision >= 128
      - id: l_late_journal_len_minus_one
        type: s4
        if: _root.file_revision >= 108
      - id: l_late_journal_len
        type: s4
        if: _root.file_revision >= 108

  journal_length_header:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_journal_data_len
        type: s4
        if: _root.file_revision >= 61

  main_compression_header_pre4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: unknown
        size: 34
        if: _root.file_revision >= 30
      - id: l_text_len
        type: s4
        if: _root.file_revision >= 30

  main_compression_header_post4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: unknown1
        size: 24
        if: _root.file_revision >= 61
      - id: l_str_len1
        type: s4
        if: _root.file_revision >= 61
      - id: l_str_len2
        type: s4
        if: _root.file_revision >= 61
      - id: unknown2
        size: 20
        if: _root.file_revision >= 61
      - id: unknown3
        size: 6
        if: _root.file_revision >= 108

  channel_compression_header:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: unknown
        size: 44
        if: _root.file_revision >= 41
      - id: l_channel_label_len
        type: s4
        if: _root.file_revision >= 41
      - id: l_unit_label_len
        type: s4
        if: _root.file_revision >= 41
      - id: l_uncompressed_len
        type: s4
        if: _root.file_revision >= 41
      - id: l_compressed_len
        type: s4
        if: _root.file_revision >= 41

  marker_header_pre4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_length
        type: s4
        if: _root.file_revision >= 30
      - id: l_markers
        type: s4
        if: _root.file_revision >= 30

  marker_header_post4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_length
        type: s4
        if: _root.file_revision >= 61
      - id: l_markers_extra
        type: s4
        if: _root.file_revision >= 61
      - id: l_markers
        type: s4
        if: _root.file_revision >= 61
      - id: unknown
        size: 6
        if: _root.file_revision >= 61
      - id: sz_defl
        size: 5
        if: _root.file_revision >= 61
      - id: unknown2
        type: s2
        if: _root.file_revision >= 61
      - id: unknown3
        size: 8
        if: _root.file_revision >= 121
      - id: unknown4
        size: 8
        if: _root.file_revision >= 128

  marker_pre_item_metadata_header_pre4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: tag
        size: 4
        if: _root.file_revision >= 30
      - id: l_item_count
        type: s4
        if: _root.file_revision >= 30
      - id: s_unknown
        size: 76
        if: _root.file_revision >= 30

  marker_item_metadata_header:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_unknown1
        type: s4
        if: _root.file_revision >= 30
      - id: l_marker_number
        type: s4
        if: _root.file_revision >= 30
      - id: b_unknown2
        size: 12
        if: _root.file_revision >= 30
      - id: rgba_color
        size: 4
        if: _root.file_revision >= 30
      - id: h_marker_tag
        type: s2
        if: _root.file_revision >= 30
      - id: h_marker_type_id
        type: s2
        if: _root.file_revision >= 30

  marker_item_header_pre4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_sample
        type: s4
        if: _root.file_revision >= 30
      - id: f_selected
        type: s2
        if: _root.file_revision >= 36
      - id: f_text_locked
        type: s2
        if: _root.file_revision >= 30
      - id: f_position_locked
        type: s2
        if: _root.file_revision >= 30
      - id: n_text_length
        type: s2
        if: _root.file_revision >= 30

  marker_item_header_post4:
    meta:
      endian:
        switch-on: _root.is_le
        cases:
          'true': le
          'false': be
    seq:
      - id: l_sample
        type: u4
        if: _root.file_revision >= 61
      - id: unknown
        size: 4
        if: _root.file_revision >= 61
      - id: n_channel
        type: s2
        if: _root.file_revision >= 61
      - id: s_marker_style
        size: 4
        if: _root.file_revision >= 61
      - id: ll_date_created
        type: u8
        if: _root.file_revision >= 128
      - id: unknown3
        size: 8
        if: _root.file_revision >= 121
      - id: n_text_length
        type: s2
        if: _root.file_revision >= 61

