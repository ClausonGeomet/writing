"""Deterministic unit checks for the additive QAR phase model."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
from contextlib import contextmanager

spec = importlib.util.spec_from_file_location('phase_models', Path(__file__).with_name('phase_models.py'))
assert spec is not None and spec.loader is not None
pm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pm)

C = np.array([[1., 0., 0.], [0., 1., 1.]])
H = np.array([[0., 1., 0.]])
CAL = np.array([[.6, .25, .15], [.65, .2, .15], [.55, .3, .15], [.58, .24, .18]])
RC = np.eye(3) * .0001


@contextmanager
def raises_value_error(message=None):
    try:
        yield
    except ValueError as exc:
        assert message is None or message in str(exc), str(exc)
    else:
        raise AssertionError('ValueError not raised')


def prior():
    return pm.fit_initial(CAL, RC)


def test_initial_covariance_and_input_immutability():
    y = CAL.copy()
    p = prior()
    assert p['status'] == 'ok'
    assert np.allclose(CAL, y)
    assert np.allclose(p['diagnostics']['raw_mean'], CAL.mean(axis=0))
    assert np.all(np.linalg.eigvalsh(p['tangent_cov']) > 0)
    assert np.allclose(np.array(p['precision']) @ np.ones(3), 0, atol=1e-12)
    assert np.isclose(sum(p['center']), 1)


def test_full_covariance_and_singular_a0_minimum_norm():
    obs = np.array([.4, .6])
    R = np.array([[.02, .006], [.006, .03]])
    a0 = pm.estimate(obs, C, R, method='A0')
    assert a0['status'] == 'ok'
    assert np.allclose(a0['mass'], [.4, .3, .3], atol=1e-10)
    assert a0['kkt_residual'] < 1e-8
    hard = pm.estimate(obs, C, R, method='A0', support=(True, False, True))
    assert np.allclose(hard['mass'], [.4, 0, .6], atol=1e-10)
    a2 = pm.estimate(np.array([.41, .59, .28]), np.vstack((C, H)),
                     np.array([[.02, .005, 0], [.005, .03, .003], [0, .003, .01]]),
                     prior(), method='A2')
    assert a2['status'] == 'ok'
    assert a2['constraint_violation'] <= 1e-9
    anatase_only = pm.estimate([.2], H, [[.01]], method='A0')
    assert np.allclose(anatase_only['mass'], [.4, .2, .4], atol=1e-10)


def test_fixed_quartz_singular_ties_and_zero_strength():
    covariance = np.array([[.0001, .00003], [.00003, .0001]])
    y = [.62, .39]
    fixed_prior = pm.fit_initial(
        [[.6, .2, .2], [.6, .3, .1], [.6, .25, .15]], RC)
    assert fixed_prior['status'] == 'ok'
    for method, kwargs in (('A0', {}), ('A2', {'prior': fixed_prior, 'strength': 0})):
        result = pm.estimate(y, C, covariance, method=method, fixed_quartz=.6, **kwargs)
        assert result['status'] == 'ok'
        assert np.allclose(result['mass'], [.6, .2, .2], atol=1e-8, rtol=0)
        assert result['constraint_violation'] <= 1e-9
        assert result['kkt_residual'] <= 1e-8
    for quartz, expected in ((0., [0., .5, .5]), (1., [1., 0., 0.])):
        result = pm.estimate(y, C, covariance, method='A0', fixed_quartz=quartz)
        assert result['status'] == 'ok'
        assert np.allclose(result['mass'], expected, atol=1e-8, rtol=0)
        assert result['kkt_residual'] <= 1e-8
    hard = pm.estimate(y, C, covariance, method='A0', fixed_quartz=.6,
                       support=(True, False, True))
    assert hard['status'] == 'ok'
    assert np.allclose(hard['mass'], [.6, 0., .4], atol=1e-8, rtol=0)


def test_fixed_quartz_full_rank_correlated_solution():
    design = np.vstack((C, H))
    covariance = np.array([[.0001, .00003, .00001],
                           [.00003, .0001, .00002],
                           [.00001, .00002, .00008]])
    y = np.array([.62, .39, .23])
    weight = np.linalg.inv(covariance)
    residual_bulk = y[:2] - np.array([.6, .4])
    expected_a = y[2] + weight[2, :2] @ residual_bulk / weight[2, 2]
    result = pm.estimate(y, design, covariance, method='A0', fixed_quartz=.6)
    assert result['status'] == 'ok'
    assert np.allclose(result['mass'], [.6, expected_a, .4 - expected_a],
                       atol=1e-8, rtol=0)
    assert result['kkt_residual'] <= 1e-8


def test_correlated_assay_closed_form_oracle():
    design = np.eye(3)
    y = np.array([.48, .37, .25])
    R = np.array([[.04, .015, 0], [.015, .03, .005], [0, .005, .02]])
    expected = y - R @ np.ones(3) * ((sum(y) - 1) / (np.ones(3) @ R @ np.ones(3)))
    result = pm.estimate(y, design, R, method='A0')
    assert result['status'] == 'ok'
    assert np.allclose(result['mass'], expected, atol=1e-10)
    assert result['kkt_residual'] <= 1e-8


def test_whole_unit_exclusion_and_cross_repeat_covariance():
    p = prior()
    stack = {'unit_id': 'one', 'D': np.vstack((C, H, H)),
             'y': np.array([.4, .6, .26, .24]), 'cov': np.eye(4) * .002}
    other = {'unit_id': 'two', 'D': C, 'y': np.array([.7, .3]), 'cov': np.eye(2) * .001}
    original = stack['y'].copy()
    all_units = pm.update_mean(p, [stack, other])
    loo = pm.update_mean(p, [stack, other], exclude_unit='one')
    only_other = pm.update_mean(p, [other])
    assert all_units['status'] == loo['status'] == 'ok'
    assert np.allclose(loo['center'], only_other['center'])
    assert np.allclose(loo['mean_cov'], only_other['mean_cov'])
    assert np.allclose(stack['y'], original)
    assert loo['diagnostics']['included_unit_ids'] == ['two']
    assert not np.allclose(all_units['mean_cov'], loo['mean_cov'])
    withheld = dict(stack, y='withheld measurement', D='withheld design', cov='withheld covariance')
    withheld_loo = pm.update_mean(p, [withheld, other], exclude_unit='one')
    assert np.allclose(withheld_loo['center'], only_other['center'])
    with raises_value_error('unknown'):
        pm.update_mean(p, [stack], exclude_unit='missing')
    with raises_value_error('duplicate'):
        pm.update_mean(p, [stack, stack])


def test_high_information_full_covariance_prior_can_be_consumed():
    initial = prior()
    design = np.vstack((C, H))
    covariance = 1e-5 * (.7 * np.eye(3) + .3 * np.ones((3, 3)))
    units = [dict(unit_id=f'unit-{i}', D=design,
                  y=np.array([.61 + .005 * i, .39 - .005 * i, .24 + .002 * i]),
                  cov=covariance) for i in range(5)]

    updated = pm.update_mean(initial, units)
    assert updated['status'] == 'ok'
    information = np.array(updated['diagnostics']['mean_information'])
    mean_cov = np.array(updated['mean_cov'])
    tangent_cov = np.array(updated['tangent_cov'])
    precision = np.array(updated['precision'])
    assert np.array_equal(information, information.T)
    assert np.array_equal(mean_cov, mean_cov.T)
    assert np.array_equal(precision, precision.T)
    assert np.allclose(information @ mean_cov, np.eye(2), atol=1e-12, rtol=0)
    expected_precision = pm.U @ np.linalg.solve(tangent_cov + mean_cov, pm.U.T)
    assert np.allclose(precision, expected_precision, atol=1e-12, rtol=0)

    center_only = pm.estimate(None, None, None, updated, method='A1')
    fitted = pm.estimate([.63, .37, .26], design, covariance, updated, method='A2')
    chained = pm.update_mean(updated, [])
    assert center_only['status'] == fitted['status'] == chained['status'] == 'ok'
    assert np.allclose(center_only['mass'], updated['center'], atol=1e-12, rtol=0)
    assert fitted['constraint_violation'] <= pm.FEAS_TOL
    assert fitted['kkt_residual'] <= pm.KKT_TOL
    assert np.array_equal(chained['diagnostics']['mean_information'], information)

    malformed_covariance = covariance.copy()
    malformed_covariance[0, 1] += 2e-12
    with raises_value_error('nonsymmetric'):
        pm.update_mean(initial, [dict(units[0], cov=malformed_covariance)])
    malformed_prior = dict(updated, diagnostics=dict(updated['diagnostics']))
    malformed_information = information.copy()
    malformed_information[0, 1] += 2e-12
    malformed_prior['diagnostics']['mean_information'] = malformed_information.tolist()
    with raises_value_error('nonsymmetric'):
        pm.estimate(None, None, None, malformed_prior, method='A1')


def test_information_and_envelope():
    R = np.eye(2) * .001
    info = pm.information(C, R)
    assert (info['rank_raw'], info['rank_augmented'], info['rank_tangent']) == (2, 2, 1)
    phase = pm.information(np.vstack((C, H)), np.eye(3) * .001)
    assert phase['rank_tangent'] == 2
    assert pm.information(C, R, support=(True, False, True))['rank_support_tangent'] == 1
    broad = pm.abundance_bounds([.4, .6], C, R)
    narrow = pm.abundance_bounds([.4, .6, .2], np.vstack((C, H)), np.eye(3) * .001)
    assert broad['status'] == narrow['status'] == 'ok'
    assert np.all(np.array(narrow['width']) <= np.array(broad['width']) + 1e-10)
    assert pm.abundance_bounds([2., 2.], C, R)['status'] == 'infeasible_envelope'


def test_validation_and_center_only():
    p = prior()
    a1 = pm.estimate(None, None, None, p, method='A1')
    assert a1['status'] == 'ok' and a1['reconstructed'] is None
    with raises_value_error():
        pm.fit_initial([[1, 0, False], [0, 1, 0]], RC)
    with raises_value_error():
        pm.estimate([.4, .6], C, np.array([[1, 2], [2, 1]]), method='A0')
    with raises_value_error():
        pm.estimate(None, None, None, p, method='A1', support=(True, False, True))
    bad = dict(p, precision=np.eye(3).tolist())
    with raises_value_error('inconsistent'):
        pm.estimate([.4, .6], C, np.eye(2), bad, method='A2')
