import { IAcmgRule } from "../../models/acmg-rule.model";
import { ISample, ISampleSummary } from "../../models/view-samples.model";

export interface ILoadAllSamplesResponse {
  isSuccess: boolean;
  sampleList?: ISampleSummary[] | null;
}

export interface IGetSampleRequest {
  sampleId: string;
}

export interface IGetSampleResponse {
  isSuccess: boolean;
  sample: ISample;
  acmgRules: IAcmgRule[];
}

export interface IAddAcmgRuleRequest {
  sampleId: string;
  variantId: number;
  ruleId: number;
}
